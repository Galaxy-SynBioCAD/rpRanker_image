import csv
import os
import itertools
import cobra
import logging
import collections
import numpy as np
import sqlite3
import pickle
import copy
import logging
import gzip

class InputReader:
    """Parser for the databases and the files related to the pathway ranking output

        WARNING: if you define inputPath, then all the files must have specific names to
        make sure that it can find the appropriate files
    """
    def __init__(self, inputPath=None):
        if inputPath and inputPath[-1:]=='/':
            inputPath = inputPath[:-1]
        #cache files
        self.globalPath = inputPath
        self.cc_preprocess = None
        self.rr_reactions = None
        self.deprecatedMNXM_mnxm = None
        self.mnxm_dG = None
        #input files
        self.rp_paths = None
        self.cofactors_rp_paths = None
        self.rp_smiles = None
        self.cobra_model = None
        self.model_chemicals = None
        self.model_compartments = None
        self.rp_transformation = None
        self.rp_smiles_inchi = None
        self.smiles_inchi = None
        if not self._loadCache(os.getcwd()+'/cache'):
            raise ValueError


    ########################### PRIVATE FUNCTIONS ###############################


    def _checkFilePath(self, path, filename):
        """Check that the directory and the filename are valid and choose to use
        either the local or the global path
        """
        if path==None:
            if self.globalPath==None:
                logging.error('Both global path and local are not set')
                return None
            else:
                if os.path.isdir(self.globalPath):
                    try:
                        fName = [i for i in os.listdir(self.globalPath) 
                                    if not i.find(filename)==-1 
                                    if not i[-3:]=='swp'
                                    if not i[-1]=='#'][0]
                        logging.info('Automatically selected '+str(fName))
                    except IndexError:
                        logging.error('Problem finding the correct '+str(filename)+' in '+str(path))
                        return None
                    if os.path.isfile(self.globalPath+'/'+fName):
                        return self.globalPath+'/'+fName
                    else:
                        logging.error('Global path file: '+str(fName)+', does not exist')
                        return None
                else:
                    logging.error('Global path is not a directory: '+str(self.globalPath))
                    return None
        else:
            if path[-1:]=='/':
                path = path[:-1]
            if os.path.isdir(path):
                if os.path.isfile(path+'/'+filename):
                    return path+'/'+filename
                else:
                    logging.error('The file is not valid: '+str(path+'/'+filename))
                    return None
            else:
                logging.error('Local path is not a directory: '+str(path))
                return None


    def _loadCache(self, path):
        """Open the cache files
        """
        try:
            self.cc_preprocess = np.load(self._checkFilePath(path, 'cc_preprocess.npz'))
        except FileNotFoundError:
            logging.error('The file '+str(path+'/cc_preprocess.npz')+' does not seem to exist')
            return False
        try:
            self.rr_reactions = pickle.load(open(self._checkFilePath(path, 'rr_reactions.pickle'), 'rb'))
        except FileNotFoundError:
            logging.error('The file '+str(path, '/rr_reactions.pickle')+' does not seem to exist')
            return False
        try:
            self.deprecatedMNXM_mnxm = pickle.load(open(self._checkFilePath(path, 'deprecatedMNXM_mnxm.pickle'), 'rb'))
        except FileNotFoundError:
            logging.error('The file '+str(path+'/deprecatedMNXM_mnxm.pickle')+' does not seem to exist')
            return False
        try:
            self.mnxm_dG = pickle.load(open(self._checkFilePath(path, 'mnxm_dG.pickle'), 'rb'))
        except FileNotFoundError:
            logging.error('The file '+str(path+'/mnxm_dG.pickle')+' does not seem to exist')
            return False
        """
        try:
            self.smiles_inchi = pickle.load(open(self._checkFilePath(path, 'smiles_inchi.pickle'), 'rb'))
        except FileNotFoundError:
            logging.error('The file '+str(path+'/smiles_inchi.pickle')+' does not seem to exists')
        """
        try:
            self.smiles_inchi = pickle.load(gzip.open(self._checkFilePath(path, 'smiles_inchi.pickle.gz'), 'rb'))
        except FileNotFoundError:
            logging.error('The file '+str(path+'/smiles_inchi.pickle')+' does not seem to exists')
        return True


    ######################## PUBLIC FUNCTIONS ################################### 


    def all(self, path=None):
        """Public function to parse all the files given a inputPath
        """
        #TODO: have the option to parse the flatfile instead of the whole database
        #TODO: have the ability to make SQL requests (such as only reactions related to particular)
        if not path==None:
            self.globalPath = path
        if self.globalPath==None:
            return False
        if not os.path.isdir(self.globalPath):
            logging.error('The global path is not a directory: '+str(self.globalPath))
            return False
        self.rp_paths = self.outPaths()
        self.rp_transformation, self.rp_smiles_inchi = self.transformation()
        self.rp_smiles = self.compounds()
        self.cofactors_rp_paths = self.addCofactors(self.rp_paths, 
                                                    self.rr_reactions, 
                                                    self.rp_smiles,
                                                    self.rp_transformation,
                                                    self.smiles_inchi,
                                                    self.rp_smiles_inchi)
        self.cobra_model = self.model()
        self.model_chemicals = self.chemicals()
        self.model_compartments = self.compartments()
        return True


    def compounds(self, path=None):
        """ Method to parse all the RP output compounds.
        TODO: remove
        """
        rp_compounds = {}
        try:
            with open(self._checkFilePath(path, 'compounds')) as f:
                reader = csv.reader(f, delimiter='\t')
                next(reader)
                for row in reader:
                    rp_compounds[row[0]] = row[1]
        except (TypeError, FileNotFoundError) as e:
            logging.error('Could not read the compounds file ('+str(path)+')')
            return {}
        return rp_compounds


    def transformation(self, path=None):
        """Extract the reaction rules from the retroPath2.0 output using the scope.csv file
        """
        rp_transformation = {}
        smiles_inchi = {}
        try:
            with open(self._checkFilePath(path, 'scope.csv')) as f:
                reader = csv.reader(f, delimiter=',')
                next(reader)
                for row in reader:
                    rp_transformation[row[1]] = row[2]
                    if not row[3] in smiles_inchi:
                        smiles_inchi[row[3]] = row[4]
                    if not row[5] in smiles_inchi:
                        smiles_inchi[row[5]] = row[6]
        except (TypeError, FileNotFoundError) as e:
            logging.error('Could not read the compounds file ('+str(path)+')')
            return {}
        return rp_transformation, smiles_inchi


    def outPaths(self, path=None):
        """RP2path Metabolic pathways form out_paths.csv
        create all the different values for heterologous paths from the RP2path out_paths.csv file
        Note that path_step are in reverse order here
        """
        ########## open either the global path or the local defined path ############
        #### (with priority with the local path)
        try:
            rp_paths = {}
            with open(self._checkFilePath(path, 'out_paths')) as f:
                reader = csv.reader(f)
                next(reader)
                current_path_id = 0
                path_step = 0
                for row in reader:
                    if not int(row[0])==current_path_id:
                        path_step = 0
                    else:
                        path_step += 1
                    current_path_id = int(row[0])
                    ################################################################
                    # WARNING: we are using ONLY the first rule and not the others #
                    ################################################################
                    #for singleRule in row[2].split(','):
                    tmpReac = {'rule_id': row[2].split(',')[0],
                            'right': {},
                            'left': {},
                            'step': path_step,
                            'path_id': int(row[0]),
                            'transformation_id': row[1][:-2]}
                    for l in row[3].split(':'):
                        tmp_l = l.split('.')
                        try:
                            #tmpReac['left'].append({'stoichio': int(tmp_l[0]), 'name': tmp_l[1]})
                            mnxm = ''
                            if tmp_l[1] in self.deprecatedMNXM_mnxm:
                                mnxm = self.deprecatedMNXM_mnxm[tmp_l[1]]
                            else:
                                mnxm = tmp_l[1]
                            tmpReac['left'][mnxm] = int(tmp_l[0])
                        except ValueError:
                            logging.error('Cannot convert tmp_l[0] to int ('+str(tmp_l[0])+')')
                            return {}
                    for r in row[4].split(':'):
                        tmp_r = r.split('.')
                        try:
                            #tmpReac['right'].append({'stoichio': int(tmp_r[0]), 'name': tmp_r[1]})
                            mnxm = ''
                            if tmp_r[1] in self.deprecatedMNXM_mnxm:
                                mnxm = self.deprecatedMNXM_mnxm[tmp_r[1]]
                            else:
                                mnxm = tmp_r[1]
                            tmpReac['right'][mnxm] = int(tmp_r[0])
                        except ValueError:
                            logging.error('Cannot convert tmp_r[0] to int ('+str(tmp_r[0])+')')
                            return {}
                    try:
                        if not int(row[0]) in rp_paths:
                            rp_paths[int(row[0])] = []
                        rp_paths[int(row[0])].insert(0, tmpReac)
                    except ValueError:
                        logging.error('Cannot convert path_id to int ('+str(row[0])+')')
                        return {}
            ####### now check where are the duplicates path_steps in each path and duplicate if yes ###
            toRet_rp_paths = [] # we make this into a list instead of a dict 
            #to find index positions in an array: usage find([1,3,4,5],[2,3])
            find = lambda searchList, elem: [[i for i, x in enumerate(searchList) if x == e] for e in elem]
            #loop through all path metabolic steps in order
            for path in rp_paths:
                dupli_items = [item for item, count in collections.Counter([i['step'] for i in rp_paths[path]]).items() if count>1]
                dupli_index = find([i['step'] for i in rp_paths[path]], dupli_items)
                flat_dupli_index = [item for sublist in dupli_index for item in sublist]
                if not dupli_items:
                    toRet_rp_paths.append(rp_paths[path])
                else:
                    keep_always_index = [i for i in [y for y in range(len(rp_paths[path]))] if i not in flat_dupli_index]
                    for dupli_include in list(itertools.product(*dupli_index)):
                        toAdd_index = list(keep_always_index+list(dupli_include))
                        new_path = []
                        for ta_i in toAdd_index:
                            new_path.append(rp_paths[path][ta_i])
                        new_path = sorted(new_path, key=lambda k: k['step'], reverse=True)
                        toRet_rp_paths.append(new_path)
            #self.rp_paths = toRet_rp_paths
            return toRet_rp_paths
        except (TypeError, FileNotFoundError) as e:
            logging.error('Could not read the out_paths file ('+str(path)+')')
            return {}


    def OLD_outPaths(self, path=None):
        """RP2path Metabolic pathways form out_paths.csv
        create all the different values for heterologous paths from the RP2path out_paths.csv file
        Note that path_step are in reverse order here
        """
        ########## open either the global path or the local defined path ############
        #### (with priority with the local path)
        try:
            rp_paths = {}
            with open(self._checkFilePath(path, 'out_paths.csv')) as f:
                reader = csv.reader(f)
                next(reader)
                current_path_id = 0
                path_step = 0
                for row in reader:
                    if not int(row[0])==current_path_id:
                        path_step = 0
                    else:
                        path_step += 1
                    current_path_id = int(row[0])
                    for singleRule in row[2].split(','):
                        tmpReac = {'rule_id': singleRule,
                                'right': {},
                                'left': {},
                                'step': path_step,
                                'path_id': int(row[0]),
                                'transformation_id': row[1][:-2]}
                        for l in row[3].split(':'):
                            tmp_l = l.split('.')
                            try:
                                #tmpReac['left'].append({'stoichio': int(tmp_l[0]), 'name': tmp_l[1]})
                                mnxm = ''
                                if tmp_l[1] in self.deprecatedMNXM_mnxm:
                                    mnxm = self.deprecatedMNXM_mnxm[tmp_l[1]]
                                else:
                                    mnxm = tmp_l[1]
                                tmpReac['left'][mnxm] = int(tmp_l[0])
                            except ValueError:
                                logging.error('Cannot convert tmp_l[0] to int ('+str(tmp_l[0])+')')
                                return {}
                        for r in row[4].split(':'):
                            tmp_r = r.split('.')
                            try:
                                #tmpReac['right'].append({'stoichio': int(tmp_r[0]), 'name': tmp_r[1]})
                                mnxm = ''
                                if tmp_r[1] in self.deprecatedMNXM_mnxm:
                                    mnxm = self.deprecatedMNXM_mnxm[tmp_r[1]]
                                else:
                                    mnxm = tmp_r[1]
                                tmpReac['right'][mnxm] = int(tmp_r[0])
                            except ValueError:
                                logging.error('Cannot convert tmp_r[0] to int ('+str(tmp_r[0])+')')
                                return {}
                        try:
                            if not int(row[0]) in rp_paths:
                                rp_paths[int(row[0])] = []
                            rp_paths[int(row[0])].insert(0, tmpReac)
                        except ValueError:
                            logging.error('Cannot convert path_id to int ('+str(row[0])+')')
                            return {}
            ####### now check where are the duplicates path_steps in each path and duplicate if yes ###
            toRet_rp_paths = [] # we make this into a list instead of a dict 
            #to find index positions in an array: usage find([1,3,4,5],[2,3])
            find = lambda searchList, elem: [[i for i, x in enumerate(searchList) if x == e] for e in elem]
            #loop through all path metabolic steps in order
            for path in rp_paths:
                dupli_items = [item for item, count in collections.Counter([i['step'] for i in rp_paths[path]]).items() if count>1]
                dupli_index = find([i['step'] for i in rp_paths[path]], dupli_items)
                flat_dupli_index = [item for sublist in dupli_index for item in sublist]
                if not dupli_items:
                    toRet_rp_paths.append(rp_paths[path])
                else:
                    keep_always_index = [i for i in [y for y in range(len(rp_paths[path]))] if i not in flat_dupli_index]
                    for dupli_include in list(itertools.product(*dupli_index)):
                        toAdd_index = list(keep_always_index+list(dupli_include))
                        new_path = []
                        for ta_i in toAdd_index:
                            new_path.append(rp_paths[path][ta_i])
                        new_path = sorted(new_path, key=lambda k: k['step'], reverse=True)
                        toRet_rp_paths.append(new_path)
            #self.rp_paths = toRet_rp_paths
            return toRet_rp_paths
        except (TypeError, FileNotFoundError) as e:
            logging.error('Could not read the out_paths file ('+str(path)+')')
            return {}


    def addCofactors(self, in_rp_paths, rr_reactions, rp_smiles, rp_transformation, smiles_inchi, rp_smiles_inchi):
        """Adds the cofactors to the retropath reactions
        """
        rp_paths = copy.deepcopy(in_rp_paths)
        for path in rp_paths:
            step_iden_CMP = {}
            for step in path:
                try:
                    reac = copy.deepcopy(rr_reactions[step['rule_id']])
                except KeyError:
                    logging.error('Cannot find rule_id: '+step['rule_id'])
                    break
                rp_path_left = None
                rp_path_right = None
                #determine where the substrate is located in the reaction of origin to determine the directionality of the reaction
                if step['rule_id'].split('_')[1] in reac['right']:
                    rp_path_right = 'right'
                    rp_path_left = 'left'
                if step['rule_id'].split('_')[1] in reac['left']:
                    rp_path_right = 'left'
                    rp_path_left = 'right'
                if step['rule_id'].split('_')[1] in reac['right'] and step['rule_id'].split('_')[1] in reac['left']:
                    logging.error('The rule_id substrate is on the left and the right of the original reaction')
                    break
                if rp_path_left==None and rp_path_right==None:
                    logging.error('Cannot find the direction from the rule compared to rule_id')
                    break
                #take the *unique* (TODO make sure) right hand side substrate to the rule
                step_iden_CMP[[i for i in step['right'].keys()][0]] = step['rule_id'].split('_')[1]
                #find the compounds to be added to the rp_paths
                toAdd_left = reac[rp_path_left]
                toAdd_right = reac[rp_path_right]
                #remove CMP from the toAdds
                for i in step_iden_CMP:
                    #WARNING: There is a small chance that you will remove elements with MNX codes that 
                    #are in fact not described as the CMP
                    toAdd_left.pop(step_iden_CMP[i], None)
                    toAdd_right.pop(step_iden_CMP[i], None)
                #add the cofactors
                for toAdd in toAdd_left:
                    if not toAdd in step['left']:
                        step['left'][toAdd] = toAdd_left[toAdd]
                for toAdd in toAdd_right:
                    if not toAdd in step['right']:
                        step['right'][toAdd] = toAdd_right[toAdd]
                #this is to add the mnx identifiers to the reactions - should not since the RP cofactors 
                #are not the same as the ones from MNX. Determine if that is valid
                #step['cmp_mnx'] = step_iden_CMP
        ##### change the format of rp_paths from lists to dictionnaries #####
        out_rp_paths = {}
        for path in rp_paths:
            out_rp_paths[path[0]['path_id']] = {}
            out_rp_paths[path[0]['path_id']]['path'] = {}
            out_rp_paths[path[0]['path_id']]['dG'] = None
            out_rp_paths[path[0]['path_id']]['dG_uncertainty'] = None
            out_rp_paths[path[0]['path_id']]['flux_sink'] = None
            out_rp_paths[path[0]['path_id']]['flux_biomass'] = None
            out_rp_paths[path[0]['path_id']]['flux_target'] = None
            out_rp_paths[path[0]['path_id']]['flux_splitObj'] = None
            out_rp_paths[path[0]['path_id']]['flux_biLevel'] = None
            out_rp_paths[path[0]['path_id']]['number_of_interventions'] = None
            for step in path:
                out_rp_paths[step['path_id']]['path'][step['step']] = {}
                out_rp_paths[step['path_id']]['path'][step['step']]['steps'] = {}
                out_rp_paths[step['path_id']]['path'][step['step']]['dG'] = None
                out_rp_paths[step['path_id']]['path'][step['step']]['dG_uncertainty'] = None
                out_rp_paths[step['path_id']]['path'][step['step']]['origin_reaction'] = step['rule_id'].split('_')[0]
                out_rp_paths[step['path_id']]['path'][step['step']]['origin_substrate'] = step['rule_id'].split('_')[1]
                out_rp_paths[step['path_id']]['path'][step['step']]['flux_biomass'] = None
                out_rp_paths[step['path_id']]['path'][step['step']]['flux_target'] = None
                out_rp_paths[step['path_id']]['path'][step['step']]['flux_splitObj'] = None
                out_rp_paths[step['path_id']]['path'][step['step']]['flux_biLevel'] = None
                try: 
                    out_rp_paths[step['path_id']]['path'][step['step']]['reaction_smiles'] = rp_transformation[step['transformation_id']]
                except KeyError:
                    out_rp_paths[step['path_id']]['path'][step['step']]['reaction_smiles'] = None
                for compound in step['left']:
                    out_rp_paths[step['path_id']]['path'][step['step']]['steps'][compound] = {}
                    out_rp_paths[step['path_id']]['path'][step['step']]['steps'][compound]['stoichiometry'] = -step['left'][compound]
                    out_rp_paths[step['path_id']]['path'][step['step']]['steps'][compound]['dG'] = None
                    out_rp_paths[step['path_id']]['path'][step['step']]['steps'][compound]['dG_uncertainty'] = None
                    ####### SMILES ######
                    try:
                        out_rp_paths[step['path_id']]['path'][step['step']]['steps'][compound]['smiles'] = rp_smiles[compound]
                    except KeyError:
                        try:
                            out_rp_paths[step['path_id']]['path'][step['step']]['steps'][compound]['smiles'] = smiles_inchi[compound]['smiles']
                        except KeyError:
                            out_rp_paths[step['path_id']]['path'][step['step']]['steps'][compound]['smiles'] = None 
                    ####### InChI #######
                    try:
                        out_rp_paths[step['path_id']]['path'][step['step']]['steps'][compound]['inchi'] = smiles_inchi[rp_smiles[compound]]
                    except KeyError:
                        try:
                            out_rp_paths[step['path_id']]['path'][step['step']]['steps'][compound]['inchi'] = smiles_inchi[compound]['inchi']
                        except KeyError:
                            out_rp_paths[step['path_id']]['path'][step['step']]['steps'][compound]['inchi'] = None
                for compound in step['right']:
                    out_rp_paths[step['path_id']]['path'][step['step']]['steps'][compound] = {}
                    out_rp_paths[step['path_id']]['path'][step['step']]['steps'][compound]['stoichiometry'] = step['right'][compound]
                    out_rp_paths[step['path_id']]['path'][step['step']]['steps'][compound]['dG'] = None
                    out_rp_paths[step['path_id']]['path'][step['step']]['steps'][compound]['dG_uncertainty'] = None
                    ######### SMILES #####
                    try: 
                        out_rp_paths[step['path_id']]['path'][step['step']]['steps'][compound]['smiles'] = rp_smiles[compound]
                    except KeyError:
                        try:
                            out_rp_paths[step['path_id']]['path'][step['step']]['steps'][compound]['smiles'] = smiles_inchi[compound]['smiles']
                        except KeyError:
                            out_rp_paths[step['path_id']]['path'][step['step']]['steps'][compound]['smiles'] = None
                    ######## InChI #######
                    try:
                        out_rp_paths[step['path_id']]['path'][step['step']]['steps'][compound]['inchi'] = rp_smiles_inchi[rp_smiles[compound]]
                    except KeyError:
                        try:
                            out_rp_paths[step['path_id']]['path'][step['step']]['steps'][compound]['inchi'] = smiles_inchi[compound]['inchi']
                        except KeyError:
                            out_rp_paths[step['path_id']]['path'][step['step']]['steps'][compound]['inchi'] = None
        return out_rp_paths


    def compartments(self, path=None):
        """ Open the compartments,tsv file from Metanetx that gives the common name to the MNX ID's
            TODO: make this optional
        """
        model_compartments = {}
        try:
            with open(self._checkFilePath(path, 'compartments')) as f:
                reader = csv.reader(f, delimiter='\t')
                next(reader)
                for row in reader:
                    model_compartments[row[0]] = {'full_name': row[1], 'short_name': row[2]}
        except (TypeError, FileNotFoundError) as e:
            return {}
        return model_compartments


    def chemicals(self, path=None):
        """ Open the chemicals.tsv file from MetaNetX that describes the sink from a model with InChI
            TODO: Replace this with a method that scans an SBML document and extracts all the chemical
            species
        """
        model_chemicals = {}
        ################# open the chemicals file #############
        try:
            with open(self._checkFilePath(path, 'chemicals')) as f:
                reader = csv.reader(f, delimiter='\t')
                for row in reader:
                    ######### chemical formula #############
                    chem_formula = None
                    if not row[3]=='':
                        chem_formula = row[3]
                    ########## mass #####################
                    mass = None
                    if not row[4]=='':
                        try:
                            mass = float(row[4])
                        except ValueError:
                            logging.error('Could not convert the mass to float ('+str(row[4])+')')
                    ########## charge ##################
                    charge = None
                    if not row[5]=='':
                        try:
                            charge = int(row[5])
                        except ValueError:
                            logging.error('Could not convert charge to int ('+str(row[5])+')')
                    ######### xref #####################
                    xref = {} #construct xref dict
                    for i in list(set([i.split(':')[0] for i in row[6].split(';')])): #unique xref db names
                        xref[i] = []
                    for i in [i.split(':') for i in row[6].split(';')]:
                        if len(i)==2:
                            xref[i[0]].append(i[1])
                    model_chemicals[row[0]] = {'name': row[1], 
                            'names': row[2].split(';'),
                            'chem_formula': chem_formula,
                            'mass': mass,
                            'charge': charge,
                            'xref': xref}
        except (TypeError, FileNotFoundError) as e:
            return {}
        return model_chemicals


    #Given the path, open the model (NOTE: only MNX models for now)
    def model(self, path=None):
        try:
            return cobra.io.read_sbml_model(self._checkFilePath(path, 'model'))
        except AttributeError:
            return None


import csv
import os
import itertools
import logging
import collections
import numpy as np
import pickle
import copy
import logging
import gzip
import sys
from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi, MolToInchiKey, AddHs

import rpSBML

## @package InputReader
#
# Documentation for the input files reader of rpFBA


## \brief Class to read all the input files
#
# Contains all the functions that read the cache files and input files to reconstruct the heterologous pathways
# To include in the input directory the following files are required:
# - chemicals.csv (MetaNetX)
# - compartments.csv (MetaNetX)
# - compounds.txt (RP2paths output)
# - out_paths.csv (RP2paths output)
# - scope.csv (RetroPath2 output)
# - xref.csv (User, if needed)
# - sink.csv (User, if needed)
class InputReader:
    """ WARNING: if you define inputPath, then all the files must have specific names to
        make sure that it can find the appropriate files
    """
    ## InputReader constructor
    # 
    #  @param self The object pointer
    #  @param inputPath The path to the folder that contains all the input/output files required
    #  @param Database The database name of the user's xref
    def __init__(self, inputPath=None, Database=None):
        if inputPath and inputPath[-1:]=='/':
            inputPath = inputPath[:-1]
        #cache files
        self.globalPath = inputPath
        self.cc_preprocess = None
        self.rr_reactions = None
        self.full_reactions = None
        self.deprecatedMNXM_mnxm = None
        #self.mnxm_dG = None
        #input files
        self.database = Database
        self.rp_paths = None
        self.cofactors_rp_paths = None
        self.rp_smiles = None
        self.cobra_model = None
        self.model_chemicals = None
        self.model_compartments = None
        self.rp_transformation = None
        self.rp_smiles_inchi = None
        self.smiles_inchi = None
        self.in_xref = None
        self.in_inchi = None
        self.convertid_inchi = {}
        self.convertid_xref = {}
        self.pub_mnx_chem_xref = None
        self.mnx_pub_chem_xref = None
        if not self._loadCache(os.getcwd()+'/cache'):
            raise ValueError


    #######################################################################
    ############################# PRIVATE FUNCTIONS ####################### 
    #######################################################################


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
            self.full_reactions = pickle.load(open(self._checkFilePath(path, 'full_reactions.pickle'), 'rb'))
        except FileNotFoundError:
            logging.error('The file '+str(path, '/full_reactions.pickle')+' does not seem to exist')
            return False
        try:
            self.deprecatedMNXM_mnxm = pickle.load(open(self._checkFilePath(path, 'deprecatedMNXM_mnxm.pickle'), 'rb'))
        except FileNotFoundError:
            logging.error('The file '+str(path+'/deprecatedMNXM_mnxm.pickle')+' does not seem to exist')
            return False
        '''
        try:
            self.mnxm_dG = pickle.load(open(self._checkFilePath(path, 'mnxm_dG.pickle'), 'rb'))
        except FileNotFoundError:
            logging.error('The file '+str(path+'/mnxm_dG.pickle')+' does not seem to exist')
            return False
        '''
        try:
            self.smiles_inchi = pickle.load(gzip.open(self._checkFilePath(path, 'smiles_inchi.pickle.gz'), 'rb'))
        except FileNotFoundError:
            logging.error('The file '+str(path+'/smiles_inchi.pickle')+' does not seem to exists')
            return False
        try:
            self.pub_mnx_chem_xref = pickle.load(gzip.open(self._checkFilePath(path, 'pub_mnx_chem_xref.gz'), 'rb'))
        except FileNotFoundError:
            logging.error('The file '+str(path+'/pub_mnx_chem_xref.pickle')+' does not seem to exists')
            return False
        try:
            self.mnx_pub_chem_xref = pickle.load(gzip.open(self._checkFilePath(path, 'mnx_pub_chem_xref.gz'), 'rb'))
        except FileNotFoundError:
            logging.error('The file '+str(path+'/mnx_pub_chem_xref.pickle')+' does not seem to exists')
            return False
        return True


    ###############################################################
    ############################# PUBLIC FUNCTIONS ################
    ############################################################### 

    ## Run all the functions
    #
    #  Run the functions required to read input and cache files to reconstruct pathways. Requires : scope.csv, chemicals.csv, ompartments.csv, compounds.txt, out_paths.csv.
    #
    #  @param self Object pointer
    #  @param path The path of input files
    def all(self, path=None):
        """Public function to parse all the files given a inputPath
        NEED TO GET RID OF THIS FUNCTION -- can provide as an example or shift to main()
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
        self.in_xref = self.user_xref()
        self.in_inchi = self.user_sink()
        self.rp_paths = self.outPaths()
        self.rp_transformation, self.rp_smiles_inchi = self.transformation()
        self.rp_smiles = self.compounds()
        #self.cobra_model = self.model()
        self.model_chemicals = self.chemicals()
        self.model_compartments = self.compartments()
        self.cofactors_rp_paths = self.addCofactors(self.rp_paths, 
                                                    self.rr_reactions,
                                                    self.full_reactions, 
                                                    self.rp_smiles,
                                                    self.rp_transformation,
                                                    self.smiles_inchi,
                                                    self.rp_smiles_inchi,
                                                    self.model_chemicals)
        return True
    

    ## Convert chemical depiction to others type of depictions
    #
    #
    #
    #  @param self The onject pointer
    #  @param idepic string depiction to be converted, str
    #  @param itype type of depiction provided as input, str
    #  @param otype types of depiction to be generated, {"", "", ..}
    #  @return odepic generated depictions, {"otype1": "odepic1", ..}
    def convert_depiction(self, idepic, itype='smiles', otype={'inchikey'}):
        """Usage example:
        - convert_depiction(idepic='CCO', otype={'inchi', 'smiles', 'inchikey'})
        - convert_depiction(idepic='InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', itype='inchi', otype={'inchi', 'smiles', 'inchikey'})
        """
        # Import (if needed)
        if itype == 'smiles':
            rdmol = MolFromSmiles(idepic, sanitize=True)
        elif itype == 'inchi':
            rdmol = MolFromInchi(idepic, sanitize=True)
        else:
            raise NotImplementedError('"{}" is not a valid input type'.format(itype))
        if rdmol is None:  # Check imprt
            raise Exception('Import error from depiction "{}" of type "{}"'.format(idepic, itype))

        # Export
        odepic = dict()
        for item in otype:
            if item == 'smiles':
                odepic[item] = MolToSmiles(rdmol)  # MolToSmiles is tricky, one mays want to check the possible options..
            elif item == 'inchi':
                odepic[item] = MolToInchi(rdmol)
            elif item == 'inchikey':
                odepic[item] = MolToInchiKey(rdmol)
            else:
                raise NotImplementedError('"{}" is not a valid output type'.format(otype))
        return odepic

    ## Function to parse a given cross references file by the user
    #
    #  Extract the identifiers for each compound
    #
    #  @param self Object pointer
    #  @param path The input file path
    #  @return xref The dictionnary of cross references
    def user_xref(self, path=None):
        """The file has to be of the form : ID   DatabaseID  Database
        """
        xref = {}
        try: 
            with open(self._checkFilePath(path, 'xref')) as f:
                reader = csv.reader(f, delimiter='\t')
                next(reader)
                for row in reader:
                    for i in range(1,len(row)-1):
                        xref[row[0]] = {'{}' .format(str(row[i+1])): row[i]}
                        i += 2
            return xref
        except(FileNotFoundError, TypeError):
            logging.warning('You did not provide an xref file, or the path is not correct')
            return None
        

    ## Function to identify the user identifier from the MNX one
    #
    #  Go through the cross references files to extract the corresponding identifiers
    #
    #  @param self Object pointer
    #  @param compound The MNX identifier to be identifiy
    #  @param db The user database containing a cross reference
    #  @return finalID 
    def cmpd_identification_xref(self, compound, db):
        tmp_id = []
        for i in self.pub_mnx_chem_xref[db]:
            if compound == self.pub_mnx_chem_xref[db][i]:
                tmp_id.append(i)
        for i in self.in_xref:
            for n in tmp_id:
                if n == self.in_xref[i][db]:
                    final_id = i
        return final_id
    
    ## Function to extract the inchi for each compound in he model of the user
    #
    #  @param self The Object pointer
    #  @return inchi_sink Dictionnary of inchis for the user's ids  
    def user_sink(self, path=None):
        try:
            with open(self._checkFilePath(path, 'sink')) as f:
                    reader = csv.reader(f, delimiter='\t')
                    next(reader)
                    inchi_sink = {}
                    for row in reader:
                        inchi_sink[row[0]]=row[1]
            return inchi_sink
        except(FileNotFoundError, TypeError):
            logging.warning('You did not provide a sink file, or the path is not correct')
            return None
        

    ## Function to identify the id of a compound using the structure comparison
    #
    #  Go through cem_prop inchis dictionnary, and the user's sink to identify the compound using inchis
    #
    #  @param self The object pointer
    #  @param compound The compound to be identified
    #  @return final_id The corresponding ID
    def cmpd_identification_inchi(self, compound):
        a = {}
        if not compound in self.convertid_inchi:
            cmpd_inchikey = self.smiles_inchi[compound]['inchikey']
            for i in self.in_inchi:
                c = self.convert_depiction(idepic=self.in_inchi[i], itype='inchi', otype={'inchikey'})
                tmp_inchikey = c['inchikey']
                if cmpd_inchikey == tmp_inchikey:
                    final_id = i
                    a[compound] = i
                elif cmpd_inchikey.split('-')[0:2] == tmp_inchikey.split('-')[0:2]:
                    logging.info('The compound '+compound+' has been identified as '+str(i)+' using the two first layers of their inchikey')
                    final_id = i
                    a[compound] = i
                elif cmpd_inchikey.split('-')[0] == tmp_inchikey.split('-')[0]:
                    logging.info('The compound '+compound+' has been identified as '+str(i)+' using the first layers of their inchikey')
                    final_id = i
                    a[compound] = i
            self.convertid_inchi.update(a)
        else:
            final_id = self.convertid_inchi[compound]
        return final_id
        
    ## Function to parse the compounds.txt file
    #
    #  Extract the smile and the structure of each compounds of RP2Path output
    #
    #  @param self Object pointer
    #  @param path The compounds.txt file path
    #  @return rp_compounds Dictionnary of smile and structure for each compound
    def compounds(self, path=None):
        """ Method to parse all the RP output compounds.
        """
        rp_compounds = {}
        try:
            with open(self._checkFilePath(path, 'compounds')) as f:
                reader = csv.reader(f, delimiter='\t')
                next(reader)
                for row in reader:
                    rp_compounds[row[0]] = {'smiles': row[1]}  #, 'structure':row[1].replace('[','').replace(']','')
        except (TypeError, FileNotFoundError) as e:
            logging.error('Could not read the compounds file ('+str(path)+')')
            return {}
        return rp_compounds

    ## Function to parse the scope.csv file
    #
    #  Extract the reaction rules from the retroPath2.0 output using the scope.csv file
    #
    #  @param self Object pointer
    #  @param path The scope.csv file path
    #  @return rp_transformation Dictionnary discribing each transformation
    #  @return smiles_inchi dictionnary describing the inchi for each smile
    def transformation(self, path=None):
        rp_transformation = {}
        smiles_inchi = {}
        try:
            with open(self._checkFilePath(path, 'results.csv')) as f:
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
    
    ## Function to parse the out_paths.csv file
    #
    #  Reading the RP2path output and extract all the information for each pathway
    #
    #  @param self Object pointer
    #  @param path The out_path.csv file path
    #  @return toRet_rp_paths Pathway object
    def outPaths(self, path=None):
        """RP2path Metabolic pathways from out_paths.csv
        create all the different values for heterologous paths from the RP2path out_paths.csv file
        Note that path_step are in reverse order here
        """
        ########## open either the global path or the local defined path ############
        #### (with priority with the local path)
        try:
            rp_paths = {}
            reactions = self.rr_reactions
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
                                mnxm = self.deprecatedMNXM_mnxm[tmp_r[1]]  #+':'+self.rr_reactions[tmpReac['rule_id']]['left']
                            else:
                                mnxm = tmp_r[1]  #+':'+self.rr_reactions[tmpReac['rule_id']]['left']
                            tmpReac['right'][mnxm] = int(tmp_r[0])
                        except ValueError:
                            logging.error('Cannot convert tmp_r[0] to int ('+str(tmp_r[0])+')')
                            return {}
                    '''##associate an MNX id for compounds annotated as CMPD ##Don't work when more than 1 compound in the left side
                    to_remove = []
                    for i in self.rr_reactions[tmpReac['rule_id']]['right'].split('.'):   
                        if not i in tmpReac['left']:
                            for n in list(tmpReac['left']):
                                if n[:4] == 'CMPD':
                                    tmpReac['left'][n+':'+i] = tmpReac['left'][n]
                                    to_remove.append(n)
                    for i in to_remove:
                        del tmpReac['left'][i]'''
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

    ## Function to reconstruct the heterologous pathway
    #
    #  Read each pathway information and RetroRules information to construct heterologous pathways and add the cofactors
    #
    #  @param self Object pointer
    #  @param in_rp_paths All pathways
    #  @param rr_reactions Dictionnary of reactions ID (rules_rall)
    #  @param full_reactions Dictionnary of reactions of origin (rxn_recipes)
    #  @param rp_smiles Dictionnary of smile and structure for each compound (compounds.txt)
    #  @param rp_transformation Dictionnary discribing each transformation
    #  @param smiles_inchi dictionnary describing the inchi for each smile (scope)
    #  @param rp_smiles_inchi dictionnary describing the inchi for each smile
    #  @return out_rp_paths Complete heterologous pathway object
    def addCofactors(self, in_rp_paths, rr_reactions, full_reactions, rp_smiles, rp_transformation, smiles_inchi, rp_smiles_inchi, model_chemicals):
        """Adds the cofactors to the retropath reactions
        """
        comp_reac_smiles = {}
        rp_paths = copy.deepcopy(in_rp_paths)
        for path in rp_paths:
            for step in path:
                try:
                    reac = copy.deepcopy(full_reactions[rr_reactions[step['rule_id']]['reaction']] )
                except KeyError:
                    logging.error('Cannot find rule_id: '+step['rule_id'])
                    sys.exit('Cannot find the rule_id. You are probably not using the last configuration of retro_rules')
                    break
                rp_path_left = None
                rp_path_right = None

                #use the relative direction in retro_rules to determine the direction of the origin reaction
                if rr_reactions[step['rule_id']]['rel_direction'] == '-1':
                    rp_path_right = 'right'
                    rp_path_left = 'left'
                if rr_reactions[step['rule_id']]['rel_direction'] == '1':
                    rp_path_right = 'left'
                    rp_path_left = 'right'
                if rp_path_left==None and rp_path_right==None:
                    logging.error('Cannot find the direction from the rule compared to rule_id')
                    break

                ##remove the main right and main left compounds from the full reaction
                toremove_l = []
                toremove_r = []
                for i in (rr_reactions[step['rule_id']][rp_path_right]).split('.'):
                    if rr_reactions[step['rule_id']]['rel_direction'] == '-1':
                        for n in reac[rp_path_left].keys():
                            if i == n and not i in toremove_l:
                                toremove_l.append(i)
                    elif rr_reactions[step['rule_id']]['rel_direction'] == '1':
                        for n in reac[rp_path_right].keys():
                            if i == n and not i in toremove_r:
                                toremove_r.append(i)
                for i in toremove_l:
                    del reac[rp_path_left][i]
                for i in toremove_r:
                    del reac[rp_path_right][i]

                toremove_l = []
                toremove_r = []
                for i in (rr_reactions[step['rule_id']][rp_path_left]).split('.'):
                    if rr_reactions[step['rule_id']]['rel_direction'] == '-1':
                        for n in reac[rp_path_right].keys():
                            if i == n and not i in toremove_r:
                                toremove_r.append(i)
                    elif rr_reactions[step['rule_id']]['rel_direction'] == '1':
                        for n in reac[rp_path_left].keys():
                            if i == n and not i in toremove_l:
                                toremove_l.append(i)
                for i in toremove_l:
                    del reac[rp_path_left][i]
                for i in toremove_r:
                    del reac[rp_path_right][i]

                #find the compounds to be added to the rp_paths according to the direction
                toAdd_left = reac[rp_path_left]
                toAdd_right = reac[rp_path_right]

                #add the cofactors according to the direction
                for toAdd in toAdd_left:
                    while toAdd in self.deprecatedMNXM_mnxm:
                        toAdd = self.deprecatedMNXM_mnxm[toAdd]
                    if not toAdd in step['left']:
                        if ('MNXM' in toAdd[:4]) and (self.in_xref is not None or self.in_inchi is not None):
                            try:
                                tmp_id = self.cmpd_identification_xref(toAdd, self.database)
                                new_toAdd = tmp_id
                            except(KeyError, UnboundLocalError, TypeError):
                                try:
                                    tmp_id = self.cmpd_identification_inchi(toAdd)
                                    new_toAdd = tmp_id
                                except(KeyError, UnboundLocalError):
                                    logging.warning("No identifier have been found for the compound "+str(toAdd))
                                except TypeError:
                                    logging.warning("You did not provide a sink file")
                        step['left'][new_toAdd+':'+toAdd] = toAdd_left[toAdd]
                for toAdd in toAdd_right:
                    while toAdd in self.deprecatedMNXM_mnxm:
                        toAdd = self.deprecatedMNXM_mnxm[toAdd]
                    if not toAdd in step['right']:
                        if ('MNXM' in toAdd[:4]) and (self.in_xref is not None or self.in_inchi is not None):
                            try:
                                tmp_id = self.cmpd_identification_xref(toAdd, self.database)
                                new_toAdd = tmp_id
                            except(KeyError, UnboundLocalError, TypeError):
                                try:
                                    tmp_id = self.cmpd_identification_inchi(toAdd)
                                    new_toAdd = tmp_id
                                except(KeyError, UnboundLocalError):
                                    logging.warning("No identifier have been found for the compound "+str(toAdd))
                                except TypeError:
                                    logging.warning("You did not provide a sink file")
                        step['right'][new_toAdd+':'+toAdd] = toAdd_right[toAdd]
                #reconstruct the complete reaction smiles by adding the smiles of cofactors
                if not step['transformation_id'] in comp_reac_smiles:
                    try:
                        scope_smiles = rp_transformation[step['transformation_id']].split('>>')
                        tmp_smiles_l=''
                        tmp_smiles_r=''
                        if rr_reactions[step['rule_id']]['rel_direction'] == '-1':
                            smiles_l = scope_smiles[0]
                            smiles_r = scope_smiles[1]
                            adding_left = toAdd_right
                            adding_right = toAdd_left
                        if rr_reactions[step['rule_id']]['rel_direction'] == '1':
                            smiles_l = scope_smiles[1]
                            smiles_r = scope_smiles[0]
                            adding_left = toAdd_left
                            adding_right = toAdd_right
                        for toAdd in adding_left:
                            try:
                                tmp = self.convert_depiction(idepic=rp_smiles_inchi[rp_smiles[toAdd]['smiles']], itype='inchi', otype={'smiles'})
                                i = 0
                                while i < int(adding_left[toAdd]):
                                    tmp_smiles_l += tmp['smiles']+'.'
                                    i += 1
                            except KeyError:
                                try:
                                    tmp = self.convert_depiction(idepic=smiles_inchi[toAdd]['inchi'], itype='inchi', otype={'smiles'})
                                    i = 0
                                    while i < int(adding_left[toAdd]):
                                        tmp_smiles_l += tmp['smiles']+'.'
                                        i += 1
                                except KeyError:
                                    tmp_smiles_l=''
                        tmp_smiles_l += smiles_l
                        for toAdd in adding_right:
                            try:
                                tmp = self.convert_depiction(idepic=rp_smiles_inchi[rp_smiles[toAdd]['smiles']], itype='inchi', otype={'smiles'})
                                i = 0
                                while i < int(adding_right[toAdd]):
                                    tmp_smiles_r += tmp['smiles']+'.'
                                    i += 1
                            except KeyError:
                                try:
                                    tmp = self.convert_depiction(idepic=smiles_inchi[toAdd]['inchi'], itype='inchi', otype={'smiles'})
                                    i = 0
                                    while i < int(adding_right[toAdd]):
                                        tmp_smiles_r += tmp['smiles']+'.'
                                        i += 1
                                except KeyError:
                                    tmp_smiles_r=''
                        tmp_smiles_r += smiles_r
                        if rr_reactions[step['rule_id']]['rel_direction'] == '-1':
                            comp_reac_smiles[step['transformation_id']] = tmp_smiles_l+str('>>')+tmp_smiles_r
                        elif rr_reactions[step['rule_id']]['rel_direction'] == '1':
                            comp_reac_smiles[step['transformation_id']] = tmp_smiles_r+str('>>')+tmp_smiles_l
                    except KeyError:
                        pass

        ######### BILAL

        sbml_paths = {} #create libsbml model for each path
        compartment_id = 'c' #TODO: define this parameter as an entry to the model
        compartment_name = 'cytoplasm' #TODO: same as above
        for path in rp_paths:
            steps = [i for i in path]
            path_id = steps[0]['path_id']
            logging.info('############ Generating an sbml object for path '+str(path_id)+' ###########')
            ###### create a libSBML model ####
            rpsbml = rpSBML.rpSBML()
            #1) create a generic Model, ie the structure and unit definitions that we will use the most
            rpsbml.genericModel('RetroPath_Pathway_'+str(path_id), 'RP_model'+str(path_id))
            '''
            model = rpsbml.createModel('RetroPath2.0 Heterologous Pathway', 'rpModel_'+str(path[0]['path_id']))
            # NOTE: the unit definitions are hard to define as input.... perhaps have a default for the moment
            # mmol_per_gDW_per_hr
            unitDef = self.createUnitDefinition(model, 'mmol_per_gDW_per_hr')
            moleUnit = self.createUnit(unitDef, libsbml.UNIT_KIND_MOLE, 1, -3, 1)
            gramUnit = self.createUnit(unitDef, libsbml.UNIT_KIND_GRAM, 1, 0, 1)
            secondUnit = self.createUnit(unitDef, libsbml.UNIT_KIND_SECOND, 1, 0, 3600)
            # kj_per_mol
            gibbsDef = self.createUnitDefinition(model, 'kj_per_mol')
            kjUnit = self.createUnit(gibbsDef, libsbml.UNIT_KIND_JOULE, 1, 3, 1)
            moleUnit = self.createUnit(gibbsDef, libsbml.UNIT_KIND_MOLE, 1, 1, 1)
            #compartment
            compartment = self.createCompartment(model, 1, compartment_name, compartment_id)
            upInfParam = self.createParameter(model, 'B_INF', float('inf'), 'kj_per_mol')
            lowInfParam = self.createParameter(model, 'B__INF', float('-inf'), 'kj_per_mol')
            '''
            ##################################
            #2) create the pathway (groups)
            rpsbml.createPathway('hetero_pathway')
            #3) find all the unique species and add them to the model
            all_meta = set([i for step in steps for lr in ['left', 'right'] for i in step[lr]])
            print(all_meta)
            for meta in list(all_meta):
                for meta2 in list(all_meta):
                    if meta.split(':')[0] == meta2.split(':')[0] and not ':' in meta and ':' in meta2:
                        print(meta, meta2)
                        all_meta.remove(meta)
                print(meta)
                cmpd_meta = meta.split(':')
                meta = cmpd_meta[0]
                #BILAL this is for you to complete, add SMILES, Inchi, etc....
                ### INCHI ### 
                try:
                    inchi = rp_smiles_inchi[rp_smiles[meta]['smiles']]
                except KeyError:
                    try:
                        inchi = smiles_inchi[meta]['inchi']
                    except KeyError:
                        inchi = None
                ### SMILES ###
                try:
                    smiles = rp_smiles[meta]['smiles']
                except KeyError:
                    try:
                        smiles = smiles_inchi[meta]['smiles']
                    except KeyError:
                        smiles = None
                ### FORMULA ###
                try:
                	formula = smiles_inchi[meta]['formula']
                except KeyError:
                	formula = ''
                ### CHARGE ###
                try:
                	charge = model_chemicals[meta]['charge']
                except KeyError:
                	charge = 0
                ### COMPOUND IDENTIFICATION ###
                #if 'TARGET' in meta[:6]:  ##'CMPD' in meta[:4] or 'TARGET' in meta[:6]
                #    tmp_meta = cmpd_meta[1]
                #else:
                #    tmp_meta = cmpd_meta[0]
                
                try:
                    tmp_xref = {}
                    if cmpd_meta[1]:
                        for db in self.mnx_pub_chem_xref[cmpd_meta[1]]:
                            tmp_xref[db] = self.mnx_pub_chem_xref[cmpd_meta[1]][db]
                    rpsbml.createSpecies(meta, tmp_xref, None, inchi, smiles, compartment_id, charge, formula)   ### [TODO] add the charge, the chemical formula and the compartment
                except (KeyError, IndexError):
                    rpsbml.createSpecies(meta, {}, None, inchi, smiles, compartment_id, charge, formula)    
            #4) add the complete reactions and their annotations
            for step in path:
                #BILAL this is for you to complete
                try:
                    reac_smiles = comp_reac_smiles[step['transformation_id']] ##rp_transformation[step['transformation_id']]
                except KeyError:
                    reac_smiles = None
                #print(step)
                rpsbml.createReaction('RetroPath_Reaction_'+str(step['step']),
                        'RP'+str(step['step']),
                        'B_INF', #only for genericModel
                        'B__INF', #only for genericModel
                        step,
                        reac_smiles,
                        compartment='c')
            #5) Optional?? Add the flux objectives. Could be in another place, TBD
            rpsbml.createFluxObj('rpFBA_obj', 'RP0', 1, True)
            sbml_paths['RP_model_'+str(path_id)] = {'model': rpsbml.document, 'flux_biomass': None, 'flux_target': None, 'flux_splitObj': None, 'flux_biLevel': None}
            #6) Test writing
            rpsbml._writeSBML('rpPath_'+str(path_id), '/home/bshahin/workspace/sbml_models/')  ## TODO change the output path to be a parameter
        return sbml_paths

    ## Function to parse the compartments.csv file
    #
    #  Parse the compartments.csv file to extract the full name and short name of the different compartments
    #
    #  @param self Object pointer 
    #  @param path The compartments.csv file path
    #  @return model_compartments Dictionnary of compartments names
    def compartments(self, path=None):
            """ Open the compartments.tsv file from Metanetx that gives the common name to the MNX ID's
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

    ## Function to parse chemicals.csv file
    #
    #  Extract different information about components
    #
    #  @param self Object pointer
    #  @param The chemicals.csv file path
    #  @return model_chemicals Dictionnary of component information
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


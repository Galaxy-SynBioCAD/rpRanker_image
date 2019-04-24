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
    def __init__(self, inputPath=None, outputPath=None, Database=None):
        if inputPath and inputPath[-1:]=='/':
            inputPath = inputPath[:-1]
            if outputPath and outputPath[-1:]=='/':
                outputPath = outputPath[:-1]
        self.globalPath = inputPath
        self.outputPath = outputPath
        #cache files
        self.cc_preprocess = None
        self.rr_reactions = None
        self.full_reactions = None
        self.deprecatedMNXM_mnxm = None
        self.mnxm_dG = None
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
        try:
            self.mnxm_dG = pickle.load(open(self._checkFilePath(path, 'mnxm_dG.pickle'), 'rb'))
        except FileNotFoundError:
            logging.error('The file '+str(path+'/mnxm_dG.pickle')+' does not seem to exist')
            return False
        try:
            self.smiles_inchi = pickle.load(gzip.open(self._checkFilePath(path, 'smiles_inchi.pickle.gz'), 'rb'))
        except FileNotFoundError:
            logging.error('The file '+str(path+'/smiles_inchi.pickle')+' does not seem to exists')
            return False
        try:
        	self.chem_xref = pickle.load(gzip.open(self._checkFilePath(path, 'chemXref.pickle.gz'), 'rb'))
        except FileNotFoundError:
            logging.error('The file '+str(path+'/chemXref.pickle.gz')+' does not seem to exists')
            return False
        try:
        	self.reac_xref = pickle.load(gzip.open(self._checkFilePath(path, 'reacXref.pickle.gz'), 'rb'))
        except FileNotFoundError:
            logging.error('The file '+str(path+'/reacXref.pickle.gz')+' does not seem to exists')
            return False
        try:
        	self.comp_xref = pickle.load(gzip.open(self._checkFilePath(path, 'compXref.pickle.gz'), 'rb'))
        except FileNotFoundError:
            logging.error('The file '+str(path+'/compXref.pickle.gz')+' does not seem to exists')
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

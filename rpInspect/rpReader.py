import csv
import os
import itertools
import logging
import collections
import pickle
import logging
import gzip
import sys
import random
from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi, MolToInchiKey, AddHs

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
class rpReader:
    """ WARNING: if you define inputPath, then all the files must have specific names to
        make sure that it can find the appropriate files
    """
    ## InputReader constructor
    # 
    #  @param self The object pointer
    #  @param inputPath The path to the folder that contains all the input/output files required
    #  @param Database The database name of the user's xref
    def __init__(self):
        #cache files
        self.deprecatedMNXM_mnxm = None
        #input files
        self.rp_paths = None
        self.rp_strc = None #These are the structures contained within the 
        #TODO: not the best strategy to have this set in this fashion -- change it
        #self.model_chemicals = None
        #self.model_compartments = None
        #TODO: not the best strategy to have this set in this fashion -- change it
        self.rp_transformation = None
        #deprecated self.rp_strc_inchi = {} #if you need to add one more then add it to mnxm_strc directly -> is that a good idea?
        self.mnxm_strc = None
        #self.in_xref = None #DEPRECATED -- to redo in a different manner:
        #self.in_inchi = None #DEPRECATED -- parameter for a function that does not exist
        #DEPRECATED -- the logic behind this is to identify user id's for compounds using
        #either a file input by the user with the user id and public id or InChI
        #Inchi is converted to its key before being compared
        self.rr_reactions = None
        if not self._loadCache():
            raise ValueError


    #######################################################################
    ############################# PRIVATE FUNCTIONS ####################### 
    #######################################################################


    ## Private function to load the required cache parameters
    #
    #
    def _loadCache(self):
        dirname = os.path.dirname(os.path.abspath( __file__ ))
        try:
            self.deprecatedMNXM_mnxm = pickle.load(open(dirname+'/cache/deprecatedMNXM_mnxm.pickle', 'rb'))
        except FileNotFoundError as e:
            logging.error(e)
            return False
        try:
            self.mnxm_strc = pickle.load(gzip.open(dirname+'/cache/mnxm_strc.pickle.gz', 'rb'))
        except FileNotFoundError as e:
            logging.error(e)
            return False
        try:
            self.rr_reactions = pickle.load(open(dirname+'/cache/rr_reactions.pickle', 'rb'))
        except FileNotFoundError as e:
            logging.error(e)
            return False
        return True


    ## Convert chemical depiction to others type of depictions
    #
    # Usage example:
    # - convert_depiction(idepic='CCO', otype={'inchi', 'smiles', 'inchikey'})
    # - convert_depiction(idepic='InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', itype='inchi', otype={'inchi', 'smiles', 'inchikey'})
    #
    #  @param self The onject pointer
    #  @param idepic string depiction to be converted, str
    #  @param itype type of depiction provided as input, str
    #  @param otype types of depiction to be generated, {"", "", ..}
    #  @return odepic generated depictions, {"otype1": "odepic1", ..}
    def _convert_depiction(self, idepic, itype='smiles', otype={'inchikey'}):
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

    ###############################################################
    ############################# PUBLIC FUNCTIONS ################
    ############################################################### 


    ''' DEPRECATED - For the moment, need to find a better way to define user xref
    ## Function to parse a given cross references file by the user
    #
    #  Extract the identifiers for each compound
    #
    #  @param self Object pointer
    #  @param path The input file path
    #  @return xref The dictionnary of cross references
    def user_xref(self, path):
        """The file has to be of the form : ID   DatabaseID  Database
        """
        self.in_xref = {}
        try:
            with open(path) as f:
                reader = csv.reader(f, delimiter='\t')
                next(reader)
                for row in reader:
                    for i in range(1,len(row)-1):
                        self.in_xref[row[0]] = {'{}' .format(str(row[i+1])): row[i]}
                        i += 2
        except(FileNotFoundError, TypeError):
            logging.warning('You did not provide an xref file, or the path is not correct')
    '''

    ''' DEPRECATED - Not sure what this function does as it is not used anywhere else
    ## Function to extract the inchi for each compound in the model of the user
    #
    #  @param self The Object pointer
    #  @return inchi_sink Dictionnary of inchis for the user's ids  
    def user_sink(self, path):
        try:
            self.in_inchi = {}
            with open(path) as f :
                reader = csv.reader(f, delimiter='\t')
                next(reader)
                for row in reader:
                    self.in_inchi[row[0]] = row[1]
        except(FileNotFoundError, TypeError):
            logging.warning('You did not provide a sink file, or the path is not correct')
    '''

    ## Function to parse the compounds.txt file
    #
    #  Extract the smile and the structure of each compounds of RP2Path output
    #  Method to parse all the RP output compounds.
    #
    #  @param self Object pointer
    #  @param path The compounds.txt file path
    #  @return rp_compounds Dictionnary of smile and structure for each compound
    def compounds(self, path):
        self.rp_strc = {}
        try:
            with open(path) as f:
                reader = csv.reader(f, delimiter='\t')
                next(reader)
                for row in reader:
                    self.rp_strc[row[0]] = {'smiles': row[1]}  #, 'structure':row[1].replace('[','').replace(']','')
                    try:
                        self.rp_strc[row[0]]['inchi'] = self.mnxm_strc[row[0]]['inchi']
                    except KeyError:
                        #try to generate them yourself by converting them directly
                        try:
                            resConv = self._convert_depiction(idepic=row[1], itype='smiles', otype={'inchi'})
                            self.rp_strc[row[0]]['inchi'] = resConv['inchi']
                        except (NotImplementedError, Exception) as e:
                            logging.warning('Could not convert the following SMILES to InChI: '+str(row[1]))
                    try:
                        self.rp_strc[row[0]]['inchikey'] = self.mnxm_strc[row[0]]['inchikey']
                        #try to generate them yourself by converting them directly
                        #TODO: consider using the inchi writing instead of the SMILES notation to find the inchikey
                    except KeyError:
                        try:
                            resConv = self._convert_depiction(idepic=row[1], itype='smiles', otype={'inchikey'})    
                            self.rp_strc[row[0]]['inchikey'] = resConv['inchikey']
                        except (NotImplementedError, Exception) as e:
                            logging.warning('Could not convert the following SMILES to InChI key: '+str(row[1]))
        except (TypeError, FileNotFoundError) as e:
            logging.error('Could not read the compounds file ('+str(path)+')')


    ## Function to parse the scope.csv file
    #
    #  Extract the reaction rules from the retroPath2.0 output using the scope.csv file
    #
    #  @param self Object pointer
    #  @param path The scope.csv file path
    #  @return rp_transformation Dictionnary discribing each transformation
    #  @return mnxm_strc dictionnary describing the inchi for each smile
    def transformation(self, path):
        try:
            with open(path) as f:
                reader = csv.reader(f, delimiter=',')
                self.rp_transformation = {}
                next(reader)
                for row in reader:
                    if not row[1] in self.rp_transformation:
                        self.rp_transformation[row[1]] = {}
                        self.rp_transformation[row[1]]['rule'] = row[2]
                        self.rp_transformation[row[1]]['ec'] = [i.replace(' ', '') for i in row[11][1:-1].split(',') if not i.replace(' ', '')=='NOEC']
                        #self.rp_transformation[row[1]]['ec'] = [i.replace(' ', '') if not i.replace(' ', '')=='NOEC' else pass for i in row[11][1:-1].split(',')]
                    '''
                    #DEPRECATED: trying to remove self.mnxm_strc for space
                    if not row[3] in self.mnxm_strc:
                        self.mnxm_strc[row[3]] = row[4]
                    if not row[5] in self.mnxm_strc:
                        self.mnxm_strc[row[5]] = row[6]
                    '''
        except FileNotFoundError:
            logging.error('Could not read the compounds file: '+str(path))


    ## Function to parse the out_paths.csv file
    #
    #  Reading the RP2path output and extract all the information for each pathway
    #  RP2path Metabolic pathways from out_paths.csv
    #  create all the different values for heterologous paths from the RP2path out_paths.csv file
    #  Note that path_step are in reverse order here
    #
    #  @param self Object pointer
    #  @param path The out_path.csv file path
    #  @maxRuleId maximal numer of rules associated with a step
    #  @return toRet_rp_paths Pathway object
    def outPaths(self, path, maxRuleIds=10):
        ########## open either the global path or the local defined path ############
        #### (with priority with the local path)
        try:
            rp_paths = {}
            #reactions = self.rr_reactions
            with open(path) as f:
                reader = csv.reader(f)
                next(reader)
                current_path_id = 0
                path_step = 1
                for row in reader:
                    try:
                        if not int(row[0])==current_path_id:
                            path_step = 1
                        else:
                            path_step += 1
                        #important to leave them in order
                        current_path_id = int(row[0])
                    except ValueError:
                        logging.error('Cannot convert path_id to int ('+str(row[0])+')')
                        return {}
                    #################################
                    ruleIds = row[2].split(',')
                    if ruleIds==None:
                        logging.error('The rulesIds is None')
                        pass
                    ###WARNING: This is the part where we select some rules over others
                    # we do it by sorting the list according to their score and taking the topx
                    if len(ruleIds)>maxRuleIds:
                        logging.warning('There are too many rules, limiting the number to random top '+str(maxRuleIds))
                        try:
                            ruleIds = [y for y,_ in sorted([(i, self.rr_reactions[i]['rule_score']) for i in ruleIds])][:maxRuleIds]
                        except KeyError:
                            logging.warning('Could not select topX due inconsistencies between rules ids and rr_reactions... selecting random instead')
                            ruleIds = random.sample(ruleIds, maxRuleIds)
                    sub_path_step = 1
                    for singleRule in ruleIds:
                        tmpReac = {'rule_id': singleRule,
                                'right': {},
                                'left': {},
                                'path_id': int(row[0]),
                                'step': path_step,
                                'sub_step': sub_path_step,
                                'transformation_id': row[1][:-2]}
                        ############ LEFT ##############
                        for l in row[3].split(':'):
                            tmp_l = l.split('.')
                            try:
                                #tmpReac['left'].append({'stoichio': int(tmp_l[0]), 'name': tmp_l[1]})
                                mnxm = '' #TODO: change this
                                if tmp_l[1] in self.deprecatedMNXM_mnxm:
                                    mnxm = self.deprecatedMNXM_mnxm[tmp_l[1]]
                                else:
                                    mnxm = tmp_l[1]
                                tmpReac['left'][mnxm] = int(tmp_l[0])
                            except ValueError:
                                logging.error('Cannot convert tmp_l[0] to int ('+str(tmp_l[0])+')')
                                return {}
                        ############## RIGHT ###########
                        for r in row[4].split(':'):
                            tmp_r = r.split('.')
                            try:
                                #tmpReac['right'].append({'stoichio': int(tmp_r[0]), 'name': tmp_r[1]})
                                mnxm = '' #TODO change this
                                if tmp_r[1] in self.deprecatedMNXM_mnxm:
                                    mnxm = self.deprecatedMNXM_mnxm[tmp_r[1]]  #+':'+self.rr_reactions[tmpReac['rule_id']]['left']
                                else:
                                    mnxm = tmp_r[1]  #+':'+self.rr_reactions[tmpReac['rule_id']]['left']
                                tmpReac['right'][mnxm] = int(tmp_r[0])
                            except ValueError:
                                logging.error('Cannot convert tmp_r[0] to int ('+str(tmp_r[0])+')')
                                return {}
                        #################################
                        if not int(row[0]) in rp_paths:
                            rp_paths[int(row[0])] = {}
                        if not int(path_step) in rp_paths[int(row[0])]:
                            rp_paths[int(row[0])][int(path_step)] = {}
                        rp_paths[int(row[0])][int(path_step)][int(sub_path_step)] = tmpReac
                        sub_path_step += 1
            self.rp_paths = rp_paths
        except (TypeError, FileNotFoundError) as e:
            logging.error(e)
            logging.error('Could not read the out_paths file ('+str(path)+') ')
            return {}



    ## Generate the sink from a given model and the 
    #
    # NOTE: this only works for MNX models, since we are parsing the id
    # TODO: change this to read the annotations and extract the MNX id's
    #
    def genSink(self, rpsbml, file_out, compartment_id='MNXC3'):
        ### open the cache ###
        cytoplasm_species = []
        for i in rpsbml.model.getListOfSpecies():
            if i.getCompartment()==compartment_id:
                cytoplasm_species.append(i)
        with open(file_out, mode='w') as f:
            writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
            writer.writerow(['Name','InChI'])
            for i in cytoplasm_species:
                res = rpsbml.readAnnotation(i.getAnnotation())
                #extract the MNX id's
                try:
                    mnx = res['metanetx.chemical'][0]
                except KeyError:
                    continue
                #mnx = i.getId().split('__')[0]
                try:
                    inchi = self.mnxm_strc[mnx]['inchi']
                except KeyError:
                    inchi = None
                if mnx and inchi:
                    writer.writerow([mnx,inchi])


if __name__ == "__main__":
    #READ THE INPUT FILES AND PASS THEM TO rpFBA
    rpreader = rpFBA.rpReader()
    rpreader.compounds(params.rp2paths_compounds)
    rpreader.transformation(params.rp2paths_scope)
    rpreader.outPaths(params.rp2paths_outPaths)
    rpcofactors = rpFBA.rpCofactors(rpreader)
    rpcofactors.pathsToSBML()
    #WRITE THE TAR.XZ
    with tarfile.open('testFBAout.tar.xz', 'w:xz') as tf:
        for rpsbml_name in rpsbml_paths:
            data = libsbml.writeSBMLToString(rpsbml_paths[rpsbml_name].document).encode('utf-8')
            fiOut = BytesIO(data)
            info = tarfile.TarInfo(rpsbml_name)
            info.size = len(data)
            tf.addfile(tarinfo=info, fileobj=fiOut)


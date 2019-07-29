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
import json
from .rpSBML import rpSBML


#######################################################
################### USER DEFINED ERROR ################
#######################################################


## Error function for the convertion of structures
#
#
class Error(Exception):
    pass


## Error function for the convertion of structures
#
#
class DepictionError(Error):
    def __init__(self, message):
        #self.expression = expression
        self.message = message



## @package InputReader
#
# Documentation for the input files reader of rpFBA

## \brief Class to read all the input files
#
# Contains all the functions that read the cache files and input files to reconstruct the heterologous pathways
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
        self.rpsbml_paths = {} #keep all the generated sbml's in this parameter
        #input files
        #TODO: remove all the rp parameters since these should not be used, 
        self.rp_strc = None #These are the structures contained within the output of rp2paths
        self.rp_transformation = None
        self.deprecatedMNXM_mnxm = None
        self.deprecatedMNXR_mnxr = None
        self.mnxm_strc = None #There are the structures from MNXM
        self.rr_reactions = None
        self.chemXref = None
        self.compXref = None
        self.rp_paths = None
        self.sbml_paths = None
        #self.reacXref = None #for the moment we are not using it, we are adding heterologous reactions
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
            self.deprecatedMNXR_mnxr = pickle.load(open(dirname+'/cache/deprecatedMNXR_mnxr.pickle', 'rb'))
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
        try:
            self.chemXref = pickle.load(gzip.open(dirname+'/cache/chemXref.pickle.gz', 'rb'))
        except FileNotFoundError as e:
            logging.error(e)
            return False
        try:
            self.compXref = pickle.load(gzip.open(dirname+'/cache/compXref.pickle.gz', 'rb'))
        except FileNotFoundError as e:
            logging.error(e)
            return False
        '''
        try:
            self.reacXref = pickle.load(gzip.open(dirname+'/cache/reacXref.pickle.gz', 'rb'))
        except FileNotFoundError as e:
            logging.error(e)
            return False
        '''
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
    ############################ RP2paths entry functions #########
    ############################################################### 

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
                        except DepictionError as e:
                            logging.warning('Could not convert the following SMILES to InChI: '+str(row[1]))
                    try:
                        self.rp_strc[row[0]]['inchikey'] = self.mnxm_strc[row[0]]['inchikey']
                        #try to generate them yourself by converting them directly
                        #TODO: consider using the inchi writing instead of the SMILES notation to find the inchikey
                    except KeyError:
                        try:
                            resConv = self._convert_depiction(idepic=row[1], itype='smiles', otype={'inchikey'})    
                            self.rp_strc[row[0]]['inchikey'] = resConv['inchikey']
                        except DepictionError as e:
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
                                'rule_score': self.rr_reactions[singleRule]['rule_score'],
                                'right': {},
                                'left': {},
                                'path_id': int(row[0]),
                                'step': path_step,
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
                        #rp_paths[int(row[0])][int(path_step)] = tmpReac
                        sub_path_step += 1
            self.rp_paths = rp_paths
        except (TypeError, FileNotFoundError) as e:
            logging.error(e)
            logging.error('Could not read the out_paths file ('+str(path)+') ')
            return {}


    ## TODO: switch this from generic to defined, with defined bounds
    #
    # rp_paths structure is the following {1: {1: {1: {'rule_id': '', 'right': {}, 'left': {}, 'path_id': int, 'step': int, 'sub_step': int, 'transformation_id': ''}, ...}, ...}, ...}
    # a single step looks like this {'rule_id': 'RR-01-503dbb54cf91-49-F', 'right': {'TARGET_0000000001': 1}, 'left': {'MNXM2': 1, 'MNXM376': 1}, 'path_id': 1, 'step': 1, 'sub_step': 1, 'transformation_id': 'TRS_0_0_17'}
    #
    #TODO: remove the default MNXC3 compartment ID
    def pathsToSBML(self, pathId='rp_pathway', compartment_id='MNXC3'):
        #if the output folder does not exist then create it
        #for path in self.rp_paths:
        #sbmlThermo = rpThermo.rpThermo()
        self.sbml_paths = {}
        for pathNum in self.rp_paths:
            #first level is the list of lists of sub_steps
            #second is itertools all possible combinations using product
            altPathNum = 1
            for comb_path in list(itertools.product(*[[y for y in self.rp_paths[pathNum][i]] for i in self.rp_paths[pathNum]])):
                steps = []
                for i in range(len(comb_path)):
                    steps.append(self.rp_paths[pathNum][i+1][comb_path[i]])
                path_id = steps[0]['path_id']
                rpsbml = rpSBML('rp_'+str(path_id)+'_'+str(altPathNum))
                #1) create a generic Model, ie the structure and unit definitions that we will use the most
                ##### TODO: give the user more control over a generic model creation:
                #   -> special attention to the compartment
                rpsbml.genericModel('RetroPath_Pathway_'+str(path_id)+'_'+str(altPathNum), 'RP_model_'+str(path_id)+'_'+str(altPathNum), self.compXref[compartment_id])
                #2) create the pathway (groups)
                rpsbml.createPathway(pathId)
                #3) find all the unique species and add them to the model
                all_meta = set([i for step in steps for lr in ['left', 'right'] for i in step[lr]])
                for meta in all_meta:
                    #here we want to gather the info from rpReader's rp_strc and mnxm_strc
                    try:
                        rpsbml.createSpecies(meta, self.chemXref[meta], None, self.rp_strc[meta]['inchi'], self.rp_strc[meta]['inchikey'], self.rp_strc[meta]['smiles'], compartment_id)
                    except KeyError:    
                        try:
                            rpsbml.createSpecies(meta, {}, None, self.rp_strc[meta]['inchi'], self.rp_strc[meta]['inchikey'], self.rp_strc[meta]['smiles'], compartment_id)
                        except KeyError:
                            logging.error('Could not create the following metabolite in either rpReaders rp_strc or mnxm_strc: '+str(meta))
                #4) add the complete reactions and their annotations
                for step in steps:
                    #add the substep to the model
                    step['sub_step'] = altPathNum
                    rpsbml.createReaction('RP'+str(step['step']), # parameter 'name' of the reaction deleted : 'RetroPath_Reaction_'+str(step['step']),
                            'B_999999', #only for genericModel
                            'B_0', #only for genericModel
                            step,
                            self.rp_transformation[step['transformation_id']]['rule'],
                            {}, #with new reactions, its impossible to find the same ones
                            self.rp_transformation,
                            compartment_id)
                #5) adding the consumption of the target
                targetRule = {'rule_id': None, 'left': {[i for i in all_meta if i[:6]=='TARGET'][0]: 1}, 'right': [], 'step': None, 'sub_step': None, 'path_id': None, 'transformation_id': None, 'rule_score': None}
                rpsbml.createReaction('targetSink',
                        'B_999999',
                        'B_0',
                        targetRule,
                        None,
                        {},
                        self.rp_transformation,
                        compartment_id,
                        True)
                #6) Optional?? Add the flux objectives. Could be in another place, TBD
                #rpsbml.createFluxObj('rpFBA_obj', 'RP0', 1, True)
                rpsbml.createFluxObj('rpFBA_obj', 'targetSink', 1, True)
                #self.sbml_paths['rp_'+str(path_id)] = rpsbml
                #self.sbml_paths['rp_'+str(step['path_id'])+'_'+str(step['sub_step'])] = rpsbml
                self.sbml_paths['rp_'+str(step['path_id'])+'_'+str(altPathNum)] = rpsbml
                altPathNum += 1


    #######################################################################
    ############################# JSON input ############################## 
    #######################################################################


    ## Function to generate an SBLM model from a json file
    #
    #  Read the json files of a folder describing pathways and generate an SBML file for each
    #
    #  @param self Object pointer
    #  @return rpsbml.document the SBML document
    #  TODO: remove the default MNXC3 compartment ID
    #  TODO: change the ID of all species to take a normal string and not sepcial caracters
    def jsonToSBML(self, jsonfile, pathId='rp_pathway', compartment_id='MNXC3'):
        pathNum = 1
        #TODO: ask how are the JSON is organised. i.e. is there one file per pathway or multiple ones
        with open(jsonfile) as json_data:
            ## Load json data as a dictionnary
            data = json.load(json_data)
            ## 1) create a generic Model, i.e the structure and unit definitions that we will be used the most
            rpsbml = rpSBML.rpSBML('rp_'+str(pathNum))
            #rpsbml.genericModel('RetroPath_Pathway_'+str(pathNum), 'RP_model'+str(pathNum)) # TODO: add self.compXref as per Melchior
            rpsbml.genericModel('RetroPath_Pathway_'+str(pathNum), 'RP_model'+str(pathNum), self.compXref[compartment_id])
            ## 2) create the pathway (groups)
            rpsbml.createPathway(pathId) ## Create the heterologous pathway
            self.rp_paths = {}
            rp_stpechio = {}
            for node in data['elements']['nodes']:
                ## ListOfSpecies
                if node['data']['type'] == 'compound':
                    try:
                        tmp_inchi = self.convert_depiction(node['data']['SMILES'], 'smiles', {'inchi'})
                    except DepictionError as e:
                        logging.warning('Cannot convert the smiles to inchi')
                        logging.warning(e)
                    ### rp_strc like (compounds.txt)
                    self.rp_strc[node['data']['id']] = {'smiles': node['data']['SMILEs'], 'inchi': tmp_inchi, 'inchikey': node['data']['id']}
                    ### create a specie
                    rpsbml.createSpecies(node['data']['id'].split('-')[0], self.inchikey_mnxm, None, tmp_inchi['inchi'], node['data']['id'], node['data']['SMILES'], compartment_id)
                    if node['data']['isSource'] == 1:
                        target_ID = node['data']['id'].split('-')[0]
                ## Create a dictionnary containing all the reactions 
                elif node['data']['type'] == 'reaction':
                    step = node['data']['iteration']
                    ### rp_transformation like (scope.scv)
                    self.rp_transformation[node['data']['id']] = {'rule': node['data']['Reaction SMILES'], 'ec': [i for i in [node['data']['EC number']]]}
                    ### rp_paths like (out_paths.csv)
                    self.rp_paths[pathNum][step][1] = {'rule_id': node['data']['Rule ID'],
                                                 'right':{},
                                                 'left':{},
                                                 'step': node['data']['iteration'], #node['data']['id'].split('-')[-1],
                                                 'EC_number': node['data']['EC number'],
                                                 'smiles': node['data']['Reaction SMILES'],
                                                 'diameter': node['data']['Diameter'],
                                                 'score': node['data']['Score'],
                                                 'iteration': node['data']['Iteration'],
                                                 'transformation_id': node['data']['id']}
                    rp_stpechio.update(node['data']['Stoechiometry'])
            ## Substrats and products for each reactions in the reactions dictionnary
            for reaction in data['elements']['edges']:
                if len(reaction['data']['target'].split('-')) == 3:
                    self.rp_paths[pathNum][reaction['data']['source'].split('-')[-1]][1]['left'][reaction['data']['target'].split('-')[0]] = rp_stpechio[reaction['data']['target']]
                else:
                    self.rp_paths[pathNum][reaction['data']['target'].split('-')[-1]][1]['right'][reaction['data']['source'].split('-')[0]] = rp_stpechio[reaction['data']['source']]
            ## ListOfReactions
            for step in self.rp_paths[pathNum]:
                rpsbml.createReaction('RetroPath_Reaction_'+step.keys().split('-')[7], # name of the reaction, to change
                        'B_999999', #only for genericModel
                        'B_0', #only for genericModel
                        step[1],
                        self.rp_transformation[step[1]['transformation_id']]['rule'],
                        {}, #self.reacXref,
                        self.rp_transformation,
                        compartment_id)
            ## targetSink reaction    
            targetRule = {'rule_id': None, 'left': {target_ID: 1}, 'right': [], 'step': None, 'path_id': None, 'transformation_id': None, 'rule_score': None}
            rpsbml.createReaction('targetSink',
                    'B_999999',
                    'B_0',
                    targetRule,
                    None,
                    {}, #
                    self.rp_transformation,
                    compartment_id,
                    True)
            rpsbml.createFluxObj('rpFBA_obj', 'targetSink', 1, True)
            self.sbml_paths['rp_'+str(pathNum)] = rpsbml
            #TODO: these will be handled by galaxy tools and not here:
            ## Writting 
            #libsbml.writeSBML(rpsbml.document, path+'/RP_model_from_json'+str(pathNum)+'.xml')
            #return rpsbml.document


    #TODO: move this to another place

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
                res = rpsbml.readMIRIAMAnnotation(i.getAnnotation())
                #extract the MNX id's
                try:
                    mnx = res['metanetx'][0]
                except KeyError:
                    continue
                #mnx = i.getId().split('__')[0]
                try:
                    inchi = self.mnxm_strc[mnx]['inchi']
                except KeyError:
                    inchi = None
                if mnx and inchi:
                    writer.writerow([mnx,inchi])


'''TODO: Need to update this test or write actual test functions
#TODO: update this thing
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
'''

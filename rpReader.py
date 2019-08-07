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
import copy

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
        #self.rpsbml_paths = {} #keep all the generated sbml's in this parameter
        #input files
        #TODO: remove all the rp parameters since these should not be used, 
        self.rp_strc = None #These are the structures contained within the output of rp2paths
        self.rp_transformation = {}
        self.deprecatedMNXM_mnxm = None
        self.deprecatedMNXR_mnxr = None
        self.mnxm_strc = None #There are the structures from MNXM
        self.inchikey_mnxm = None #There are the mnxmIDs for InChIkeys
        self.rr_reactions = None
        self.chemXref = None
        self.compXref = None
        self.rp_paths = {}
        self.sbml_paths = {}
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
            self.inchikey_mnxm = pickle.load(gzip.open(dirname+'/cache/inchikey_mnxm.pickle.gz', 'rb'))
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
                rpsbml.genericModel('RetroPath_Pathway_'+str(path_id)+'_'+str(altPathNum), 
                        'RP_model_'+str(path_id)+'_'+str(altPathNum), 
                        self.compXref[compartment_id])
                #2) create the pathway (groups)
                rpsbml.createPathway(pathId)
                #3) find all the unique species and add them to the model
                all_meta = set([i for step in steps for lr in ['left', 'right'] for i in step[lr]])
                for meta in all_meta:
                    #here we want to gather the info from rpReader's rp_strc and mnxm_strc
                    try:
                        rpsbml.createSpecies(meta, 
                                self.chemXref[meta], 
                                None, 
                                self.rp_strc[meta]['inchi'], 
                                self.rp_strc[meta]['inchikey'], 
                                self.rp_strc[meta]['smiles'], 
                                compartment_id)
                    except KeyError:    
                        try:
                            rpsbml.createSpecies(meta, 
                                    {}, 
                                    None, 
                                    self.rp_strc[meta]['inchi'], 
                                    self.rp_strc[meta]['inchikey'], 
                                    self.rp_strc[meta]['smiles'], 
                                    compartment_id)
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
                            compartment_id,
                            self.rp_transformation[step['transformation_id']]['rule'],
                            self.rp_transformation[step['transformation_id']]['ec'])
                    #5) adding the consumption of the target
                targetStep = {'rule_id': None, 'left': {[i for i in all_meta if i[:6]=='TARGET'][0]: 1}, 'right': [], 'step': None, 'sub_step': None, 'path_id': None, 'transformation_id': None, 'rule_score': None}
                rpsbml.createReaction('targetSink',
                        'B_999999',
                        'B_0',
                        targetStep,
                        compartment_id)
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


    ## Pass a dictionnary of JSON dict and convert them to SBML
    #
    #
    def collectionJSON(self, collJson): 
        pathNum = 1
        sub_path = 1
        #According to the rule applied there should be subpaths ==> in this case we ignore that
        self.sbml_paths = {}
        for i in collJson:
            rpsbml = self.jsonToSBML(collJson[i], pathNum)
            self.sbml_paths['rp_'+str(pathNum)+'_'+str(sub_path)] = rpsbml
            pathNum += 1
            

    ## Function to generate an SBLM model from a json file
    #
    #  Read the json files of a folder describing pathways and generate an SBML file for each
    #
    #  @param self Object pointer
    #  @return rpsbml.document the SBML document
    #  TODO: remove the default MNXC3 compartment ID
    #  TODO: change the ID of all species to take a normal string and not sepcial caracters
    # WARNING: We are only using a single rule (technically with the highest diameter)
    def jsonToSBML(self, json_dict, pathNum=1, pathId='rp_pathway', compartment_id='MNXC3'):
        #pathNum = 1
        #create the SBML
        rpsbml = rpSBML('rp_'+str(pathNum))
        rpsbml.genericModel('RetroPath_Pathway_'+str(pathNum), 'RP_model'+str(pathNum), self.compXref[compartment_id])
        rpsbml.createPathway(pathId) 
        ### gather the data
        all_reac = {}
        rp_paths = {}
        stochio = {}
        inchikey_cid = {}
        for node in json_dict['elements']['nodes']:
            # add the species to the SBML directly
            if node['data']['type']=='compound':
                try:
                    #hack, pick the smallest mnx ID if there are multiple ones
                    cid = sorted(self.inchikey_mnxm[node['data']['id']]['mnx'], key=lambda x: int(x[4:]))[0]
                except KeyError:
                    cid = node['data']['id'].replace('-', '') 
                inchikey_cid[node['data']['id'].replace('-', '')] = cid
                try:
                    rpsbml.createSpecies(cid, 
                                        self.chemXref[cid],
                                        None, 
                                        node['data']['InChI'], 
                                        node['data']['id'], 
                                        node['data']['SMILES'], 
                                        compartment_id)
                except KeyError:
                    rpsbml.createSpecies(cid, 
                                        {},
                                        None, 
                                        node['data']['InChI'], 
                                        node['data']['id'], 
                                        node['data']['SMILES'], 
                                        compartment_id)
            ## Create a dictionnary containing all the reactions 
            elif node['data']['type']=='reaction':
                #note sure about this... need to ask 
                step = node['data']['Iteration']
                #hack, pick the first rule, usually have the biggest diameter
                r_id = sorted(node['data']['Rule ID'], key=lambda x: int(x.split('-')[-2]), reverse=True)[0]
                all_reac[node['data']['id']] = {'rule_id': r_id, 
                        'right': {},
                        'left': {},
                        'path_id': pathNum,
                        'step': step,
                        'sub_step': 1,
                        'transformation_id': node['data']['id'],
                        'rule_score': node['data']['Score'],
                        'smiles': node['data']['Reaction SMILES'],
                        'ec': list(filter(None, [i for i in node['data']['EC number']]))}
                stochio[node['data']['id']] = node['data']['Stoechiometry']
        for reaction_node in json_dict['elements']['edges']:
            if not len(reaction_node['data']['source'].split('-'))==3:
                if not reaction_node['data']['source'] in all_reac:
                    logging.warning('The following reaction was not found in the JSON elements: '+str(reaction_node['data']['source']))
                else:
                    rid = reaction_node['data']['source']
                    try:
                        cid = inchikey_cid[reaction_node['data']['target'].replace('-', '')]
                    except KeyError:
                        cid = reaction_node['data']['target'].replace('-', '')
                    try:
                        all_reac[rid]['right'][cid] = stochio[rid][cid]
                    except KeyError:
                        all_reac[rid]['right'][cid] = 1.0
            #target
            if not len(reaction_node['data']['target'].split('-'))==3:
                if not reaction_node['data']['target'] in all_reac:
                    logging.warning('The following reaction was not found in the JSON elements: '+str(reaction_node['data']['source']))
                else:
                    rid = reaction_node['data']['target']
                    try:
                        cid = inchikey_cid[reaction_node['data']['source'].replace('-', '')]
                    except KeyError:
                        cid = reaction_node['data']['source'].replace('-', '')
                    try:
                        all_reac[rid]['left'][cid] = stochio[rid][cid]
                    except KeyError:
                        all_reac[rid]['left'][cid] = 1.0
        # now that all the information has been gathered pass it to the SBML
        for rule_id in all_reac:
            rpsbml.createReaction('RP'+str(all_reac[rule_id]['step']), 
                    'B_999999', #only for genericModel
                    'B_0', #only for genericModel
                    #{'left': all_reac[rule_id]['left'], 'right': all_reac[rule_id]['right']},
                    all_reac[rule_id],
                    compartment_id,
                    all_reac[rule_id]['smiles'],
                    all_reac[rule_id]['ec'],
                    {})
        return rpsbml

    #############################################################################################
    ############################### validation data tsv #########################################
    #############################################################################################


    ##
    #
    #
    def parseValidation(self, inFile):
        data = {}
        try:
            for row in csv.DictReader(open(inFile), delimiter='\t'):
                ######## pathId ######
                try:
                    pathID = int(row['pathway_ID'])
                except ValueError:
                    logging.error('Cannot convert pathway ID: '+str(row['pathway_ID']))
                    continue
                if not pathID in data:
                    data[pathID] = {}
                    data[pathID]['isValid'] = True
                    data[pathID]['steps'] = {}
                ####### step #########
                try:
                    stepID = int(row['step'])
                except ValueError:
                    logging.error('Cannot convert step ID: '+str(row['step']))
                    data[pathID]['isValid'] = False
                    continue
                if stepID==0:
                    continue
                elif stepID==1:
                    data[pathID]['organism'] = row['organism'].replace(' ', '')
                    data[pathID]['reference'] = row['reference'].replace(' ', '')
                data[pathID]['steps'][stepID] = {}
                ##### substrates #########
                data[pathID]['steps'][stepID]['substrates'] = []
                lenDBref = len(row['substrate_dbref'].split(';'))
                for i in row['substrate_dbref'].split(';'):
                    if i=='':
                        lenDBref -= 1
                lenStrc = len(row['substrate_structure'].split(';'))
                for i in row['substrate_structure'].split(';'):
                    if i=='':
                        lenStrc -= 1
                lenSub = len(row['substrate_name'].split(';'))
                for i in row['substrate_name'].split(';'):
                    if i=='':
                        lenSub -= 1
                if lenSub==lenStrc==lenSub:
                    for name, inchi, dbrefs in zip(row['substrate_name'].split(';'), 
                            row['substrate_structure'].split(';'), 
                            row['substrate_dbref'].split(';')):
                        tmp = {}
                        tmp['inchi'] = inchi.replace(' ', '')
                        tmp['name'] = name
                        tmp['dbref'] = {}
                        for dbref in dbrefs.split('|'):
                            if len(dbref.split(':'))==2:
                                db_name = dbref.split(':')[0].replace(' ', '').lower()
                                db_cid = dbref.split(':')[1].replace(' ', '')
                                if not db_name in tmp['dbref']:
                                    tmp['dbref'][db_name] = []
                                tmp['dbref'][db_name].append(db_cid)
                            else:
                                logging.warning('Ignoring the folowing product dbref ('+str(name)+'): '+str(dbref))
                                data[pathID]['isValid'] = False
                        data[pathID]['steps'][stepID]['substrates'].append(tmp)
                else:
                    logging.warning('Not equal length between substrate names, their structure or dbref ('+str(name)+'): '+str(row['substrate_name'])+' <--> '+str(row['substrate_structure'])+' <--> '+str(row['substrate_dbref']))
                    data[pathID]['isValid'] = False
                    continue
                ##### products #########
                data[pathID]['steps'][stepID]['products'] = []
                lenDBref = len(row['product_dbref'].split(';'))
                for i in row['product_dbref'].split(';'):
                    if i=='':
                        lenDBref -= 1
                lenStrc = len(row['product_structure'].split(';'))
                for i in row['product_structure'].split(';'):
                    if i=='':
                        lenStrc -= 1
                lenSub = len(row['product_name'].split(';'))
                for i in row['product_name'].split(';'):
                    if i=='':
                        lenSub -= 1
                if lenSub==lenStrc==lenDBref:
                    for name, inchi, dbrefs in zip(row['product_name'].split(';'), 
                            row['product_structure'].split(';'), 
                            row['product_dbref'].split(';')):
                        tmp = {}
                        tmp['inchi'] = inchi.replace(' ', '')
                        tmp['name'] = name
                        tmp['dbref'] = {}
                        for dbref in dbrefs.split('|'):
                            if len(dbref.split(':'))==2:
                                db_name = dbref.split(':')[0].replace(' ', '').lower()
                                db_cid = dbref.split(':')[1].replace(' ', '')
                                if not db_name in tmp['dbref']:
                                    tmp['dbref'][db_name] = []
                                tmp['dbref'][db_name].append(db_cid)
                            else:
                                data[pathID]['isValid'] = False
                                logging.warning('Ignoring the folowing product dbref ('+str(name)+'): '+str(dbref))
                        data[pathID]['steps'][stepID]['products'].append(tmp)
                else:
                    logging.warning('Not equal length between substrate names, their structure or dbref ('+str(name)+'): '+str(row['product_name'])+' <--> '+str(row['product_structure'])+' <--> '+str(row['product_dbref']))
                    data[pathID]['isValid'] = False
                data[pathID]['steps'][stepID]['ec_numbers'] = [i.replace(' ', '') for i in row['EC_number'].split(';')]
                data[pathID]['steps'][stepID]['enzyme_id'] = [i.replace(' ', '') for i in row['enzyme_identifier'].split(';')]
                data[pathID]['steps'][stepID]['enzyme_name'] = row['enzyme_name'].split(';')
        except FileNotFoundError:
            logging.error('Cannot open the file: '+str(inFile))
        #now loop through all of them and remove the invalid paths
        toRet = copy.deepcopy(data)
        for path_id in data.keys():
            if toRet[path_id]['isValid']==False:
                del toRet[path_id]
            else:
                del toRet[path_id]['isValid']
        return toRet


    ##
    #
    #
    def validationToSBML(self, inFile, compartment_id='MNXC3'): 
        data = self.parseValidation(inFile)
        self.sbml_paths = {}
        #TODO: need to exit at this loop
        for path_id in data:
            rpsbml = rpSBML('measured_'+str(path_id))
            #1) create a generic Model, ie the structure and unit definitions that we will use the most
            ##### TODO: give the user more control over a generic model creation:
            #   -> special attention to the compartment
            rpsbml.genericModel('measured_'+str(path_id), 'measured_'+str(path_id), self.compXref[compartment_id])
            #find all the chemical species and add them to an SBML
            #2) create the pathway (groups)
            rpsbml.createPathway(path_id)
            #3) find all the unique species and add them to the model
            allChem = []
            for stepNum in data[path_id]['steps']:
                #because of the nature of the input we need to remove duplicates
                for i in data[path_id]['steps'][stepNum]['substrates']+data[path_id]['steps'][stepNum]['products']:
                    if not i in allChem:
                        allChem.append(i)
            #add them to the SBML
            for chem in allChem:
                #PROBLEM: as it stands one expects the meta to be MNX
                if 'mnx' in chem['dbref']:
                    #must list the different models
                    meta = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
                else:
                    logging.warning('All species must be referenced by a MNX id or will be ignored')
                    break
                #try to conver the inchi into the other structures
                smiles = None
                inchikey = None
                try:
                    resConv = self._convert_depiction(idepic=chem['inchi'], itype='inchi', otype={'smiles','inchikey'})
                    smiles = resConv['smiles']
                    inchikey = resConv['inchikey']
                except DepictionError as e:
                    logging.warning('Could not convert the following SMILES to InChI: '+str(row[1]))
                #create a new species
                try:
                    rpsbml.createSpecies(meta, 
                            self.chemXref[meta], 
                            None, 
                            chem['inchi'], 
                            inchikey, 
                            smiles, 
                            compartment_id)
                except KeyError:
                    try:
                        rpsbml.createSpecies(meta, 
                                {}, 
                                None, 
                                chem['inchi'], 
                                inchikey, 
                                smiles, 
                                compartment_id)
                    except KeyError:
                        logging.error('Could not create the following metabolite: '+str(meta))
                        break
            #4) add the complete reactions and their annotations
            #need to convert the validation to step for reactions
            for stepNum in data[path_id]['steps']:
                toSend = {'left': {}, 'right': {}, 'rule_id': None, 'rule_score': None, 'path_id': path_id, 'step': stepNum, 'sub_step': None}
                for chem in data[path_id]['steps'][stepNum]['substrates']:
                    if 'mnx' in chem['dbref']:
                        meta = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
                        toSend['left'][meta] = 1
                    else:
                        logging.error('Need all the species to have a MNX ID')
                        break
                for chem in data[path_id]['steps'][stepNum]['products']:
                    if 'mnx' in chem['dbref']:
                        meta = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
                        toSend['right'][meta] = 1
                    else:
                        logging.error('Need all the species to have a MNX ID')
                        break
                #if all are full add it
                rpsbml.createReaction('M'+str(stepNum),
                        'B_999999', #only for genericModel
                        'B_0', #only for genericModel
                        toSend,
                        compartment_id,
                        None,
                        data[path_id]['steps'][stepNum]['ec_numbers'])
                rpsbml.createFluxObj('rpFBA_obj', 'M'+str(max(data[path_id]['steps'])), 1, True)
            self.sbml_paths['measured_'+str(path_id)] = rpsbml


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

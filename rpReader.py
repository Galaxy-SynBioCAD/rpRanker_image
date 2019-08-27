import csv
import os
import itertools
import collections
import pickle
import gzip
import sys
import random
from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi, MolToInchiKey, AddHs
import json
import copy

#from .setup_self.logger import self.logger
import logging
from .rpSBML import rpSBML


## @package rpReader
#
# Collection of functions that convert the outputs from various sources to the SBML format (rpSBML) for further analyses


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


#WARNING: if you define inputPath, then all the files must have specific names to
#make sure that it can find the appropriate files


## Class to read all the input files
#
# Contains all the functions that read the cache files and input files to reconstruct the heterologous pathways
class rpReader:
    ## InputReader constructor
    #
    #  @param self The object pointer
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self.logger.info('Starting instance of rpReader')
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
    #  Opens the previously generated cache to the object memory
    #
    # @param The oject pointer
    # @return Boolean detemining the success of the function or not
    def _loadCache(self):
        dirname = os.path.dirname(os.path.abspath( __file__ ))
        try:
            self.deprecatedMNXM_mnxm = pickle.load(open(dirname+'/cache/deprecatedMNXM_mnxm.pickle', 'rb'))
        except FileNotFoundError as e:
            self.logger.error(e)
            return False
        try:
            self.deprecatedMNXR_mnxr = pickle.load(open(dirname+'/cache/deprecatedMNXR_mnxr.pickle', 'rb'))
        except FileNotFoundError as e:
            self.logger.error(e)
            return False
        try:
            self.mnxm_strc = pickle.load(gzip.open(dirname+'/cache/mnxm_strc.pickle.gz', 'rb'))
        except FileNotFoundError as e:
            self.logger.error(e)
            return False
        ## could be deleted in favor of mnxm_strc
        try:
            self.inchikey_mnxm = pickle.load(gzip.open(dirname+'/cache/inchikey_mnxm.pickle.gz', 'rb'))
        except FileNotFoundError as e:
            self.logger.error(e)
            return False
        try:
            self.rr_reactions = pickle.load(open(dirname+'/cache/rr_reactions.pickle', 'rb'))
        except FileNotFoundError as e:
            self.logger.error(e)
            return False
        try:
            self.chemXref = pickle.load(gzip.open(dirname+'/cache/chemXref.pickle.gz', 'rb'))
        except FileNotFoundError as e:
            self.logger.error(e)
            return False
        try:
            self.compXref = pickle.load(gzip.open(dirname+'/cache/compXref.pickle.gz', 'rb'))
        except FileNotFoundError as e:
            self.logger.error(e)
            return False
        try:
            self.nameCompXref = pickle.load(gzip.open(dirname+'/cache/nameCompXref.pickle.gz', 'rb'))
        except FileNotFoundError as e:
            self.logger.error(e)
            return False
        ''' Not used as of now
        try:
            self.reacXref = pickle.load(gzip.open(dirname+'/cache/reacXref.pickle.gz', 'rb'))
        except FileNotFoundError as e:
            self.logger.error(e)
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
                            self.logger.warning('Could not convert the following SMILES to InChI: '+str(row[1]))
                    try:
                        self.rp_strc[row[0]]['inchikey'] = self.mnxm_strc[row[0]]['inchikey']
                        #try to generate them yourself by converting them directly
                        #TODO: consider using the inchi writing instead of the SMILES notation to find the inchikey
                    except KeyError:
                        try:
                            resConv = self._convert_depiction(idepic=row[1], itype='smiles', otype={'inchikey'})
                            self.rp_strc[row[0]]['inchikey'] = resConv['inchikey']
                        except DepictionError as e:
                            self.logger.warning('Could not convert the following SMILES to InChI key: '+str(row[1]))
        except (TypeError, FileNotFoundError) as e:
            self.logger.error('Could not read the compounds file ('+str(path)+')')
            return False


    ## Function to parse the scope.csv file
    #
    #  Extract the reaction rules from the retroPath2.0 output using the scope.csv file
    #
    #  @param self Object pointer
    #  @param path The scope.csv file path
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
            self.logger.error('Could not read the compounds file: '+str(path))
            return False


    #TODO: make sure that you account for the fact that each reaction may have multiple associated reactions
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
    #def outPaths(self, path, maxRuleIds=10):
    def outPathsToSBML(self, path, maxRuleIds=10, pathId='rp_pathway', compartment_id='MNXC3'):
        #try:
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
                    self.logger.error('Cannot convert path_id to int ('+str(row[0])+')')
                    #return {}
                    return False
                #################################
                ruleIds = row[2].split(',')
                if ruleIds==None:
                    self.logger.error('The rulesIds is None')
                    #pass # or continue
                    continue
                ###WARNING: This is the part where we select some rules over others
                # we do it by sorting the list according to their score and taking the topx
                tmp_rr_reactions = {}
                for r_id in ruleIds:
                    for rea_id in self.rr_reactions[r_id]:
                        tmp_rr_reactions[str(r_id)+'__'+str(rea_id)] = self.rr_reactions[r_id][rea_id]
                if len(ruleIds)>maxRuleIds:
                    self.logger.warning('There are too many rules, limiting the number to random top '+str(maxRuleIds))
                    try:
                        #ruleIds = [y for y,_ in sorted([(i, self.rr_reactions[i]['rule_score']) for i in ruleIds])][:maxRuleIds]
                        ruleIds = [y for y,_ in sorted([(i, tmp_rr_reactions[i]['rule_score']) for i in tmp_rr_reactions])][:maxRuleIds]
                    except KeyError:
                        self.logger.warning('Could not select topX due inconsistencies between rules ids and rr_reactions... selecting random instead')
                        ruleIds = random.sample(tmp_rr_reactions, maxRuleIds)
                else:
                    ruleIds = tmp_rr_reactions
                sub_path_step = 1
                for singleRule in ruleIds:
                    tmpReac = {'rule_id': singleRule.split('__')[0],
                            'mnxr': singleRule.split('__')[1],
                            'rule_score': self.rr_reactions[singleRule.split('__')[0]][singleRule.split('__')[1]]['rule_score'],
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
                            self.logger.error('Cannot convert tmp_l[0] to int ('+str(tmp_l[0])+')')
                            #return {}
                            return False
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
                            self.logger.error('Cannot convert tmp_r[0] to int ('+str(tmp_r[0])+')')
                            return {}
                    #################################
                    if not int(row[0]) in rp_paths:
                        rp_paths[int(row[0])] = {}
                    if not int(path_step) in rp_paths[int(row[0])]:
                        rp_paths[int(row[0])][int(path_step)] = {}
                    rp_paths[int(row[0])][int(path_step)][int(sub_path_step)] = tmpReac
                    #rp_paths[int(row[0])][int(path_step)] = tmpReac
                    sub_path_step += 1
        #except (TypeError, FileNotFoundError) as e:
        #    self.logger.error(e)
        #    self.logger.error('Could not read the out_paths file ('+str(path)+') ')
        #    return {}
        #### pathToSBML ####
        #self.rp_paths = rp_paths
        try:
            mnxc = self.nameCompXref[compartment_id]
        except KeyError:
            self.logger.error('Could not Xref compartment_id ('+str(compartment_id)+')')
            return False
        self.sbml_paths = {}
        #for pathNum in self.rp_paths:
        for pathNum in rp_paths:
            #first level is the list of lists of sub_steps
            #second is itertools all possible combinations using product
            altPathNum = 1
            #for comb_path in list(itertools.product(*[[y for y in self.rp_paths[pathNum][i]] for i in self.rp_paths[pathNum]])):
            for comb_path in list(itertools.product(*[[y for y in rp_paths[pathNum][i]] for i in rp_paths[pathNum]])):
                steps = []
                for i in range(len(comb_path)):
                    #steps.append(self.rp_paths[pathNum][i+1][comb_path[i]])
                    steps.append(rp_paths[pathNum][i+1][comb_path[i]])
                path_id = steps[0]['path_id']
                rpsbml = rpSBML('rp_'+str(path_id)+'_'+str(altPathNum))
                #1) create a generic Model, ie the structure and unit definitions that we will use the most
                ##### TODO: give the user more control over a generic model creation:
                #   -> special attention to the compartment
                rpsbml.genericModel('RetroPath_Pathway_'+str(path_id)+'_'+str(altPathNum),
                        'RP_model_'+str(path_id)+'_'+str(altPathNum),
                        self.compXref[mnxc],
                        compartment_id)
                #2) create the pathway (groups)
                rpsbml.createPathway(pathId)
                #3) find all the unique species and add them to the model
                all_meta = set([i for step in steps for lr in ['left', 'right'] for i in step[lr]])
                for meta in all_meta:
                    #here we want to gather the info from rpReader's rp_strc and mnxm_strc
                    try:
                        rpsbml.createSpecies(meta,
                                compartment_id,
                                self.chemXref[meta],
                                None,
                                self.rp_strc[meta]['inchi'],
                                self.rp_strc[meta]['inchikey'],
                                self.rp_strc[meta]['smiles'])
                    except KeyError:
                        try:
                            rpsbml.createSpecies(meta,
                                    compartment_id,
                                    {},
                                    None,
                                    self.rp_strc[meta]['inchi'],
                                    self.rp_strc[meta]['inchikey'],
                                    self.rp_strc[meta]['smiles'])
                        except KeyError:
                            self.logger.error('Could not create the following metabolite in either rpReaders rp_strc or mnxm_strc: '+str(meta))
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
                targetStep = {'rule_id': None, 
                        'left': {[i for i in all_meta if i[:6]=='TARGET'][0]: 1}, 
                        'right': [], 
                        'step': None, 
                        'sub_step': None, 
                        'path_id': None, 
                        'transformation_id': None, 
                        'rule_score': None, 
                        'mnxr': None}
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
        return True

    
    ## Function to group all the functions for parsing RP2 output to SBML files
    #
    # Takes RP2paths's compounds.txt and out_paths.csv and RetroPaths's *_scope.csv files and generates SBML
    #
    # @param compounds string path to RP2paths out_paths file
    # @param scope string path to RetroPaths2's scope file output
    # @param outPaths string path to RP2paths out_paths file
    # @param maxRuleIds int The maximal number of members in a single substep (Reaction Rule)
    # @param compartment_id string The ID of the SBML's model compartment where to add the reactions to
    # @return Boolean The success or failure of the function
    def rp2ToSBML(self, compounds, scope, outPaths, maxRuleIds=10, pathId='rp_pathway', compartment_id='MNXC3'):
        self.compounds(compounds)
        self.transformation(scope)
        return self.outPathsToSBML(outPaths, maxRuleIds, pathId, compartment_id)

    #######################################################################
    ############################# JSON input ##############################
    #######################################################################

    #WARNING: we are not setting any restrictions on the number of steps allowed here and instead we
    #take the rule with the highest dimension. Also assume that there is a single rule at a maximal
    #dimension
    ## Function to generate an SBLM model from a JSON file
    #
    #  Read the json files of a folder describing pathways and generate an SBML file for each
    #  TODO: remove the default MNXC3 compartment ID
    #  TODO: change the ID of all species to take a normal string and not sepcial caracters
    #  WARNING: We are only using a single rule (technically with the highest diameter)
    #
    #  @param self Object pointer
    # @param colJson Dictionnary of 
    #  @return rpsbml.document the SBML document
    def jsonToSBML(self, collJson, pathId='rp_pathway', compartment_id='MNXC3'):
        #global parameters used for all parameters
        pathNum = 1
        rp_paths = {}
        reac_smiles = {}
        reac_ecs = {}
        species_list = {}
        reactions_list = {}
        source_cid = {}
        source_stochio = {}
        cid_inchikey = {}
        sink_species = {}
        ############################################
        ############## gather the data #############
        ############################################
        for json_dict in collJson:
            ########### construct rp_paths ########
            reac_smiles[pathNum] = {}
            reac_ecs[pathNum] = {}
            species_list[pathNum] = {}
            reactions_list[pathNum] = {}
            cid_inchikey[pathNum] = {}
            sink_species[pathNum] = {}
            stochio = {}
            inchikey_cid = {}
            source_species = []
            skip_pathway = False
            for node in collJson[json_dict]['elements']['nodes']:
                ##### compounds ####
                if node['data']['type']=='compound':
                    cid_inchikey[pathNum][node['data']['id'].replace('-', '')] = node['data']['id']
                    species_list[pathNum][node['data']['id'].replace('-', '')] = {'inchi': node['data']['InChI'], 
                                    'inchikey': node['data']['id'],
                                    'smiles': node['data']['SMILES']}
                    if int(node['data']['isSource'])==1:
                        #TODO: there should always be only one source, to check
                        source_species.append(node['data']['id'].replace('-', '')) 
                        source_cid[pathNum] = node['data']['id'].replace('-', '')
                    if int(node['data']['inSink'])==1:
                        sink_species[pathNum][node['data']['id'].replace('-', '')] = node['data']['id']
                ###### reactions ######
                elif node['data']['type']=='reaction':
                    #NOTE: pick the rule with the highest diameter
                    r_id = sorted(node['data']['Rule ID'], key=lambda x: int(x.split('-')[-2]), reverse=True)[0]
                    reactions_list[pathNum][node['data']['id']] = {'rule_id': r_id,
                        'mnxr': None,
                        'right': {},
                        'left': {},
                        'path_id': pathNum,
                        'step': None,
                        'sub_step': None,
                        'transformation_id': None,
                        'rule_score': node['data']['Score']}
                    reac_smiles[pathNum][r_id] = node['data']['Reaction SMILES']
                    reac_ecs[pathNum][r_id] = list(filter(None, [i for i in node['data']['EC number']]))
                    stochio[node['data']['id']] = {}
                    for i in node['data']['Stoechiometry']:
                        stochio[node['data']['id']][i.replace('-', '')] = node['data']['Stoechiometry'][i]
            ############ make main pathway ###########
            main_reac = {}
            for reaction_node in collJson[json_dict]['elements']['edges']:
                if not len(reaction_node['data']['source'].split('-'))==3:
                    if not reaction_node['data']['source'] in reactions_list[pathNum]:
                        self.logger.error('The following reaction was not found in the JSON elements: '+str(reaction_node['data']['source']))
                        skip_pathway = True
                        break
                    else:
                        rid = reaction_node['data']['source']
                        try:
                            cid = inchikey_cid[reaction_node['data']['target'].replace('-', '')]
                        except KeyError:
                            cid = reaction_node['data']['target'].replace('-', '')
                        try:
                            reactions_list[pathNum][rid]['left'][cid] = stochio[rid][cid]
                        except KeyError:
                            reactions_list[pathNum][rid]['left'][cid] = 1.0
                            self.logger.warning('The cid ('+str(cid)+') has not been detected by stochio. Setting to 1.0')
                if not len(reaction_node['data']['target'].split('-'))==3:
                    if not reaction_node['data']['target'] in reactions_list[pathNum]:
                        self.logger.error('The following reaction was not found in the JSON elements: '+str(reaction_node['data']['source']))
                        skip_pathway = True
                        break
                    else:
                        rid = reaction_node['data']['target']
                        try:
                            cid = inchikey_cid[reaction_node['data']['source'].replace('-', '')]
                        except KeyError:
                            cid = reaction_node['data']['source'].replace('-', '')
                        try:
                            reactions_list[pathNum][rid]['right'][cid] = stochio[rid][cid]
                        except KeyError:
                            reactions_list[pathNum][rid]['right'][cid] = 1.0
                            self.logger.warning('The cid ('+str(cid)+') has not been detected by stochio. Setting to 1.0')
            ################# calculate the steps associated with the reactions_list ######
            #find the source in the LAST reaction in the pathway
            #NOTE: this assumes that the source is contained in a single final reaction and nowhere else
            last_rid = None
            step_num = len(reactions_list[pathNum])
            found_rid = []
            toFind_rid = list(reactions_list[pathNum].keys())
            for rid in reactions_list[pathNum]:    
                if all([True if i in source_species else False for i in reactions_list[pathNum][rid]['right']]):
                    reactions_list[pathNum][rid]['step'] = step_num
                    step_num -= 1
                    last_rid = rid
                    try:
                        source_stochio[pathNum] = stochio[rid][source_cid[pathNum]]
                    except KeyError:
                        source_stochio[pathNum] = 1.0
                    found_rid.append(rid)
                    toFind_rid.remove(rid)
                    break
            for rid in toFind_rid[:]:
                if all([True if i in reactions_list[pathNum][last_rid]['left'] else False for i in reactions_list[pathNum][rid]['right']]):
                    reactions_list[pathNum][rid]['step'] = step_num
                    step_num -= 1
                    last_rid = rid
                    found_rid.append(rid)
                    toFind_rid.remove(rid)
            if not toFind_rid==[]:
                self.logger.error('There are reactions unaccounted for: '+str(toFind_rid))
                skip_pathway = True
                break
            ############# find all the alternative reactions associated with a reaction rule ###
            if not skip_pathway:
                rp_paths[pathNum] = {}
                for rid in reactions_list[pathNum]:
                    rp_paths[pathNum][reactions_list[pathNum][rid]['step']] = {}
                    sub_step = 1
                    for reac_id in self.rr_reactions[reactions_list[pathNum][rid]['rule_id']]:
                        tmpReac = copy.deepcopy(reactions_list[pathNum][rid])
                        tmpReac['mnxr'] = reac_id
                        tmpReac['sub_step'] = sub_step
                        rp_paths[pathNum][reactions_list[pathNum][rid]['step']][sub_step] = tmpReac
                        sub_step += 1
            else:
                self.logger.warning('Skipping pathway '+str(pathNum))
            pathNum += 1
        ######################################
        ########### create the SBML's ########
        ######################################
        try:
            mnxc = self.nameCompXref[compartment_id]
        except KeyError:
            self.logger.error('Could not Xref compartment_id ('+str(compartment_id)+')')
            return False
        self.sbml_paths = {}
        for pathNum in rp_paths:
            #first level is the list of lists of sub_steps
            #second is itertools all possible combinations using product
            altPathNum = 1
            for comb_path in list(itertools.product(*[[y for y in rp_paths[pathNum][i]] for i in rp_paths[pathNum]])):
                steps = []
                for i in range(len(comb_path)):
                    steps.append(rp_paths[pathNum][i+1][comb_path[i]])
                path_id = steps[0]['path_id']
                rpsbml = rpSBML('rp_'+str(path_id)+'_'+str(altPathNum))
                #1) create a generic Model, ie the structure and unit definitions that we will use the most
                ##### TODO: give the user more control over a generic model creation:
                #   -> special attention to the compartment
                rpsbml.genericModel('RetroPath_Pathway_'+str(path_id)+'_'+str(altPathNum),
                        'RP_model_'+str(path_id)+'_'+str(altPathNum),
                        self.compXref[mnxc],
                        compartment_id)
                #2) create the pathway (groups)
                rpsbml.createPathway(pathId)
                #3) find all the unique species and add them to the model
                meta_to_cid = {}
                for meta in species_list[pathNum]:
                    #### beofre adding it to the model check to see if you can recover some MNXM from inchikey
                    #NOTE: only for the sink species do we try to convert to MNXM
                    if meta in list(sink_species[pathNum].keys()):
                        try:
                            #take the smallest MNX, usually the best TODO: review this
                            cid = sorted(self.inchikey_mnxm[sink_species[pathNum][meta]], key=lambda x: int(x[4:]))[0]
                            meta_to_cid[meta] = cid
                        except KeyError:
                            logging.error('Cannot find sink compound: '+str(meta))
                            return False
                    else:
                        cid = meta
                    #here we want to gather the info from rpReader's rp_strc and mnxm_strc
                    try:
                        rpsbml.createSpecies(cid,
                                compartment_id,
                                self.chemXref[cid],
                                None,
                                species_list[pathNum][meta]['inchi'],
                                species_list[pathNum][meta]['inchikey'],
                                species_list[pathNum][meta]['smiles'])
                    except KeyError:
                        try:
                            rpsbml.createSpecies(meta,
                                    compartment_id,
                                    {},
                                    None,
                                    species_list[pathNum][meta]['inchi'],
                                    species_list[pathNum][meta]['inchikey'],
                                    species_list[pathNum][meta]['smiles'])
                        except KeyError:
                            self.logger.error('Could not create the following metabolite in either rpReaders rp_strc or mnxm_strc: '+str(meta))
                #4) add the complete reactions and their annotations
                for step in steps:
                    #add the substep to the model
                    step['sub_step'] = altPathNum
                    #### try to replace the sink compounds inchikeys with mnxm
                    for direc in ['right', 'left']:
                        step_mnxm = {}
                        for meta in step[direc]:
                            try:
                                step_mnxm[meta_to_cid[meta]] = step[direc][meta]
                            except KeyError:
                                step_mnxm[meta] = step[direc][meta] 
                        step[direc] = step_mnxm
                    rpsbml.createReaction('RP'+str(step['step']),
                            'B_999999', #only for genericModel
                            'B_0', #only for genericModel
                            step,
                            compartment_id,
                            reac_smiles[pathNum][step['rule_id']],
                            reac_ecs[pathNum][step['rule_id']])
                #5) adding the consumption of the target
                targetStep = {'rule_id': None, 
                        'left': {source_cid[pathNum]: source_stochio[pathNum]}, 
                        'right': {}, 
                        'step': None, 
                        'sub_step': None, 
                        'path_id': None, 
                        'transformation_id': None, 
                        'rule_score': None, 
                        'mnxr': None}
                rpsbml.createReaction('targetSink',
                        'B_999999',
                        'B_0',
                        targetStep,
                        compartment_id)
                #6) Optional?? Add the flux objectives. Could be in another place, TBD
                rpsbml.createFluxObj('rpFBA_obj', 'targetSink', 1, True)
                self.sbml_paths['rp_'+str(step['path_id'])+'_'+str(altPathNum)] = rpsbml
                altPathNum += 1


    #############################################################################################
    ############################### validation data tsv #########################################
    #############################################################################################


    ## Function to parse the TSV of measured heterologous pathways to SBML
    #
    # Given the TSV of measured pathways, parse them to a dictionnary, readable to next be parsed
    # to SBML
    #
    # @param self object pointer
    # @param inFile The input JSON file
    # @param mnxHeader Reorganise the results around the target MNX products
    # @return Dictionnary of SBML
    def parseValidation(self, inFile, mnxHeader=False):
        data = {}
        try:
            for row in csv.DictReader(open(inFile), delimiter='\t'):
                ######## pathId ######
                try:
                    pathID = int(row['pathway_ID'])
                except ValueError:
                    self.logger.error('Cannot convert pathway ID: '+str(row['pathway_ID']))
                    continue
                if not pathID in data:
                    data[pathID] = {}
                    data[pathID]['isValid'] = True
                    data[pathID]['steps'] = {}
                ####### step #########
                try:
                    stepID = int(row['step'])
                except ValueError:
                    self.logger.error('Cannot convert step ID: '+str(row['step']))
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
                                self.logger.warning('Ignoring the folowing product dbref ('+str(name)+'): '+str(dbref))
                                data[pathID]['isValid'] = False
                        data[pathID]['steps'][stepID]['substrates'].append(tmp)
                else:
                    self.logger.warning('Not equal length between substrate names, their structure or dbref ('+str(name)+'): '+str(row['substrate_name'])+' <--> '+str(row['substrate_structure'])+' <--> '+str(row['substrate_dbref']))
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
                                self.logger.warning('Ignoring the folowing product dbref ('+str(name)+'): '+str(dbref))
                        data[pathID]['steps'][stepID]['products'].append(tmp)
                else:
                    self.logger.warning('Not equal length between substrate names, their structure or dbref ('+str(name)+'): '+str(row['product_name'])+' <--> '+str(row['product_structure'])+' <--> '+str(row['product_dbref']))
                    data[pathID]['isValid'] = False
                data[pathID]['steps'][stepID]['ec_numbers'] = [i.replace(' ', '') for i in row['EC_number'].split(';')]
                data[pathID]['steps'][stepID]['enzyme_id'] = [i.replace(' ', '') for i in row['enzyme_identifier'].split(';')]
                data[pathID]['steps'][stepID]['enzyme_name'] = row['enzyme_name'].split(';')
        except FileNotFoundError:
            self.logger.error('Cannot open the file: '+str(inFile))
        #now loop through all of them and remove the invalid paths
        toRet = copy.deepcopy(data)
        for path_id in data.keys():
            if toRet[path_id]['isValid']==False:
                del toRet[path_id]
            else:
                del toRet[path_id]['isValid']
        #reorganise the results around the target products mnx
        if not mnxHeader:
            return toRet
        else:
            toRetTwo = {}
            for pathId in toRet:
                final_pro_mnx = toRet[pathId]['steps'][max(toRet[pathId]['steps'])]['products'][0]['dbref']['mnx'][0]
                if not final_pro_mnx in toRetTwo:
                    toRetTwo[final_pro_mnx] = {}
                toRetTwo[final_pro_mnx][pathId] = toRet[pathId]
            return toRetTwo


    ## Parse the validation TSV to SBML
    #
    # Parse the TSV file to SBML format and adds them to the self.sbml_paths
    #
    # @param self Object pointer
    # @param inFile Input file
    # @param compartment_id compartment of the
    def validationToSBML(self, inFile, compartment_id='MNXC3'):
        data = self.parseValidation(inFile)
        self.sbml_paths = {}
        #TODO: need to exit at this loop
        for path_id in data:
            try:
                mnxc = self.nameCompXref[compartment_id]
            except KeyError:
                self.logger.error('Could not Xref compartment_id ('+str(compartment_id)+')')
                return False
            rpsbml = rpSBML('measured_'+str(path_id))
            #1) create a generic Model, ie the structure and unit definitions that we will use the most
            ##### TODO: give the user more control over a generic model creation:
            #   -> special attention to the compartment
            rpsbml.genericModel('measured_'+str(path_id),
                                'measured_'+str(path_id), 
                                self.compXref[mnxc], 
                                compartment_id)
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
                    #TODO: add the species with other types of xref in annotation
                    self.logger.warning('Some species are not referenced by a MNX id or will be ignored')
                    #break
                    
                #try to conver the inchi into the other structures
                smiles = None
                inchikey = None
                try:
                    resConv = self._convert_depiction(idepic=chem['inchi'], itype='inchi', otype={'smiles','inchikey'})
                    smiles = resConv['smiles']
                    inchikey = resConv['inchikey']
                except DepictionError as e:
                    self.logger.warning('Could not convert the following SMILES to InChI: '+str(row[1]))
                #create a new species
                try:
                    rpsbml.createSpecies(meta,
                            compartment_id,
                            self.chemXref[meta],
                            None,
                            chem['inchi'],
                            inchikey,
                            smiles)
                except KeyError:
                    try:
                        rpsbml.createSpecies(meta,
                                compartment_id,
                                {},
                                None,
                                chem['inchi'],
                                inchikey,
                                smiles)
                    except KeyError:
                        self.logger.error('Could not create the following metabolite: '+str(meta))
                        break
            #4) add the complete reactions and their annotations
            #need to convert the validation to step for reactions
            for stepNum in data[path_id]['steps']:
                toSend = {'left': {}, 'right': {}, 'rule_id': None, 'mnxr': None, 'rule_score': None, 'path_id': path_id, 'step': stepNum, 'sub_step': None}
                for chem in data[path_id]['steps'][stepNum]['substrates']:
                    if 'mnx' in chem['dbref']:
                        meta = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
                        toSend['left'][meta] = 1
                    else:
                        self.logger.error('Need all the species to have a MNX ID')
                        break
                for chem in data[path_id]['steps'][stepNum]['products']:
                    if 'mnx' in chem['dbref']:
                        meta = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
                        toSend['right'][meta] = 1
                    else:
                        self.logger.error('Need all the species to have a MNX ID')
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

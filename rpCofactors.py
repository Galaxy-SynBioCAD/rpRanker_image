import pickle
import os
import sys
import gzip
import copy
import logging
import itertools

from .rpSBML import rpSBML
#from rpThermo import rpThermo

#TODO: inherit InputReader instead of passing it as a parameter

## Class to add the cofactors to a monocomponent reaction to construct the full reaction
#
#
class rpCofactors:
    ## Init method
    # Here we want to seperate what is the use input and what is parsed by the cache to make sure that
    # everything is not hadled by a single 
    #
    # @param rpReader input reader object with the parsed user input and cache files required
    #DEPRECATED def __init__(self, rpReader, userXrefDbName=None):
    def __init__(self, rpsbml):
        #TODO: check that the sbml is valid
        self.rpsbml = rpsbml
        ##### stuff to load from cache #####
        self.full_reactions = None
        self.rr_reactions = None
        self.chemXref = None
        self.reacXref = None
        self.compXref = None
        if not self._loadCache():
            raise ValueError


    ######################################################
    ################## PRIVATE FUNCTIONS #################
    ######################################################


    ## Load the cache required for this class
    #
    #
    def _loadCache(self):
        dirname = os.path.dirname(os.path.abspath( __file__ ))
        try:
            self.full_reactions = pickle.load(open(dirname+'/cache/full_reactions.pickle', 'rb'))
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
        try:
            self.reacXref = pickle.load(gzip.open(dirname+'/cache/reacXref.pickle.gz', 'rb'))
        except FileNotFoundError as e:
            logging.error(e)
            return False
        return True

    
    ################################################################
    ######################### PUBLIC FUNCTIONS #####################
    ################################################################


    ## Function that converts an rpSBML model to the dictionnary that is used by the algorithm
    #
    #
    def convertSBMLdict():
        

    ## Function that converts a dictionnary used by the algorithm to rpSBML model
    #
    def convertDictSBML():


    ## Given a dictionnary describing a monocomponent reaction, add the cofactors by comparing it with the original reaction
    #
    # @param step Dictionnary describing the reaction
    # @param reac_side String 'right' or 'left' describing the direction of the monocomponent reaction compared with the original reaction
    # @param rr_reac Dictionnary describing the monocomponent reaction from RetroRules
    # @param f_reac Dictionnary describing the full original reaction
    # @param pathway_cmp_mnxm Dictionnary used to retreive the public ID of the intermediate compounds. Resets for each individual pathway
    #
    def completeReac(self, step, reac_side, rr_reac, f_reac, pathway_cmp_mnxm):
        ######## COFACTORS #####
        if reac_side=='right':
            #from the monocomponent side, remove the main species from RR
            noMain_fullReac = {i:f_reac[i] for i in f_reac if i!=self.rr_reactions[step['rule_id']]['subs_id']}
            #find the intermediate compound
            cmp_diff = list(step[reac_side].keys()-noMain_fullReac.keys())
            #add it to the dictionnary
            if len(cmp_diff)==1:
                pathway_cmp_mnxm.update({cmp_diff[0]: self.rr_reactions[step['rule_id']]['subs_id']})
            else:
                logging.warning('There are more than 1 or None of cmp_diff: '+str(cmp_diff))
        elif reac_side=='left':
            #identtify the main compounds and remove them from the full reaction
            try:
                toRem = [pathway_cmp_mnxm[i] for i in list(step[reac_side].keys()-f_reac.keys())]
                noMain_fullReac = {i:f_reac[i] for i in f_reac if i not in toRem}
            except KeyError:
                logging.warning('could not find internediate compound name')
                raise KeyError
        else:
            logging.warning('Direction can only be right or left')
            raise KeyError
        #calculate the difference between the two
        diff = {i: noMain_fullReac[i] for i in noMain_fullReac.keys()-step[reac_side].keys()}
        #update the reaction
        step[reac_side].update(diff)
        ########## STOCHIO #####
        #TODO: in both stochio and reaction rule reconstruction, if an error occurs we do not remove the step from the final
        #results ==>  consider raising the error if that is the case, depends on how important it is
        #update the stochiometry from the main reaction, including the intermediate compounds
        #pathway_mnxm_cmp = {v: k for k, v in pathway_cmp_mnxm.items()}
        for i in step[reac_side]:
            try:
                diff_stochio = f_reac[i]-step[reac_side][i]
                if diff_stochio>0:
                    if i in diff:
                        diff[i] += diff_stochio
                    else:
                        diff[i] = diff_stochio
                step[reac_side][i] = f_reac[i]
            except KeyError:
                try:
                    diff_stochio = f_reac[pathway_cmp_mnxm[i]]-step[reac_side][i]
                    if diff_stochio>0:
                        if i in diff:
                            diff[i] += diff_stochio
                        else:
                            diff[i] = diff_stochio
                    step[reac_side][i] = f_reac[pathway_cmp_mnxm[i]]
                except KeyError:
                    logging.warning('Could not find the intermediate compound in full reaction: '+str(i))
                    pass
        ######### REACTION RULE ########
        #take all the added chemical compounds, return their SMILES and add them to the appropriate side
        for i in diff:
            #for y in range(step[reac_side][i]): #based on the stochio, we would want to add as many as the original reaction specifies
            for y in range(diff[i]): #based on the stochio, we would want to add as many as the original reaction specifies
                if i in self.rpReader.mnxm_strc:
                    if 'smiles' in self.rpReader.mnxm_strc[i] and not self.rpReader.mnxm_strc[i]['smiles']==None:
                        reac_smiles = step['reaction_rule'].split('>>')
                        if reac_side=='left':
                            reac_smiles[1] += '.'+self.rpReader.mnxm_strc[i]['smiles']
                        else: #if it is anything else than left here, above should detect error
                            reac_smiles[0] += '.'+self.rpReader.mnxm_strc[i]['smiles']
                        step['reaction_rule'] = reac_smiles[0]+'>>'+reac_smiles[1]
                    else:
                        logging.warning('There are no SMILES defined for '+str(i)+' in self.rpReader.mnxm_strc[i]')
                        continue
                else:
                    logging.warning('Cannot find '+str(i)+' in self.rpReader.mnxm_strc')
                    continue


    ## Add the cofactors to monocomponent reactions
    #
    # @param step Step in a pathway
    # @param pathway_cmp_mnxm Dictionnary of intermediate compounds with their public ID's
    # @return Boolean determine if the step is to be added
    def addCofactors_step(self, step, pathway_cmp_mnxm):
        try:
            step['reaction_rule'] = self.rpReader.rp_transformation[step['transformation_id']]['rule']
            step['rule_score'] = self.rpReader.rr_reactions[step['rule_id']]['rule_score']
            if self.rr_reactions[step['rule_id']]['rel_direction']==-1:
                self.completeReac(step, 'right', self.rr_reactions[step['rule_id']]['left'], self.full_reactions[self.rr_reactions[step['rule_id']]['reac_id']]['right'], pathway_cmp_mnxm)
                self.completeReac(step, 'left', self.rr_reactions[step['rule_id']]['right'], self.full_reactions[self.rr_reactions[step['rule_id']]['reac_id']]['left'], pathway_cmp_mnxm)
            if self.rr_reactions[step['rule_id']]['rel_direction']==1:
                self.completeReac(step, 'right', self.rr_reactions[step['rule_id']]['left'], self.full_reactions[self.rr_reactions[step['rule_id']]['reac_id']]['left'], pathway_cmp_mnxm)
                self.completeReac(step, 'left', self.rr_reactions[step['rule_id']]['right'], self.full_reactions[self.rr_reactions[step['rule_id']]['reac_id']]['right'], pathway_cmp_mnxm)
        except KeyError:
            return False
        return True


    ## Function to reconstruct the heterologous pathway
    #
    #  Read each pathway information and RetroRules information to construct heterologous pathways and add the cofactors
    #
    #  @param self Object pointer
    #  @return out_rp_paths Complete heterologous pathway object
    def addCofactors(self):
        self.rp_paths = copy.deepcopy(self.rpReader.rp_paths)
        #here we return the list of keys of the pathway to avoid RuntimeError due to dictionnary size changes 
        for pathNum in list(self.rp_paths.keys()):
            toDel = []
            #'This keeps the IDs conversions to the pathway
            #We reverse the loop to ID the intermediate CMP to their original ones
            pathway_cmp_mnxm = {}
            for stepNum in sorted(list(self.rp_paths[pathNum].keys()), reverse=True):
                for sub_stepNum in sorted(list(self.rp_paths[pathNum][stepNum].keys()), reverse=True):
                    isCofacAdd = self.addCofactors_step(self.rp_paths[pathNum][stepNum][sub_stepNum], pathway_cmp_mnxm)
                    if not isCofacAdd:
                        #toDel.append((pathNum, stepNum, sub_stepNum))
                        toDel.append((stepNum, sub_stepNum))
            #delete elements that have failed in the cofactor_step method
            #reverse such that we can delete without KeyErrors
            for i in sorted(toDel, reverse=True):
                #TODO: need to do a better job at recording the pathway that is ignored here
                logging.warning('Deleting rp_paths at position: ['+str(pathNum)+']['+str(i[0])+']['+str(i[1])+']')
                try:
                    #if a step in a pathway contains only one step that failed, then delete the whole pathway
                    if len(self.rp_paths[pathNum][i[0]])==1:
                        logging.warning('There is only one instance of the subpath: ['+str(pathNum)+']['+str(i[0])+'], and it failed thus deleting the whole pathway')
                        del self.rp_paths[pathNum]
                        break
                    else:
                        #delete the sub_step
                        del self.rp_paths[pathNum][i[0]][i[1]]
                    #TODO: compare the subsequences and delete the ones that are the same -- need effecient algorithm
                except KeyError:
                    logging.warning('It seems like the position: ['+str(pathNum)+']['+str(i[0])+']['+str(i[1])+'] does not exist')
                    #TODO: better error handling here
                    return False


    ## TODO: switch this from generic to defined, with defined bounds
    #
    #
    #TODO: remove the default MNXC3 compartment ID
    def pathsToSBML(self, compartment_id='MNXC3'):
        #if the output folder does not exist then create it
        #for path in rp_paths:
        #sbmlThermo = rpThermo.rpThermo()
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
                #print(self.rp_paths)
                #print('############################')
                rpsbml.genericModel('RetroPath_Pathway_'+str(path_id)+'_'+str(altPathNum), 'RP_model_'+str(path_id)+'_'+str(altPathNum), self.compXref)
                #2) create the pathway (groups)
                rpsbml.createPathway('rp_pathway')
                #3) find all the unique species and add them to the model
                all_meta = set([i for step in steps for lr in ['left', 'right'] for i in step[lr]])
                ''' DEPRECATED -- need to find a better way to have user xref
                print(all_meta)
                for meta in list(all_meta):
                    for meta2 in list(all_meta):
                        if meta.split(':')[0] == meta2.split(':')[0] and not ':' in meta and ':' in meta2:
                            all_meta.remove(meta)
                for meta in all_meta:
                    if not ':' in meta and not 'CMPD' in meta[:4] and not 'TARGET' in meta[:6]:
                        try:
                            meta = meta+':'+self.cmpd_identification_xref(meta, self.userXrefDbName)
                        except(KeyError, UnboundLocalError):
                            logging.warning("No identifier has been found for the compound "+str(meta))
                            meta = meta
                    cmpd_meta = meta.split(':')
                    meta = cmpd_meta[0]
                '''
                for meta in all_meta:
                    #here we want to gather the info from rpReader's rp_strc and mnxm_strc
                    try:
                        rpsbml.createSpecies(meta, self.chemXref, None, self.rpReader.rp_strc[meta]['inchi'], self.rpReader.rp_strc[meta]['inchikey'], self.rpReader.rp_strc[meta]['smiles'], compartment_id)
                    except KeyError:    
                        try:
                            rpsbml.createSpecies(meta, self.chemXref, None, self.rpReader.mnxm_strc[meta]['inchi'], self.rpReader.mnxm_strc[meta]['inchikey'], self.rpReader.mnxm_strc[meta]['smiles'], compartment_id)
                        except KeyError:
                            logging.error('Could not create the following metabolite in either rpReaders rp_strc or mnxm_strc: '+str(meta))
                    ''' No need for this since at the rpCache we try to perform that conversion for MNXM, and the same is 
                    performed for intermediate compounds at the rpReader level
                    ### Attempt to generate a SMILES from Inchi and vice versa ###
                    #### and INCHI Key ###
                    if not smiles==None and inchi==None:
                        convertRes = self.convert_depiction(idepic=smiles, itype='smiles', otype={'inchi', 'inchikey'})
                        inchi = convertRes['inchi']
                        inchiKey = convertRes['inchikey']
                    if not inchi==None and smiles==None:
                        convertRes = self.convert_depiction(idepic=inchi, itype='inchi', otype={'smiles', 'inchikey'})
                        smiles = convertRes['smiles']
                        inchiKey = convertRes['inchikey']
                    '''
                    ''' DEPRECATED
                    #TODO: add at some point but as of now not very important
                    ### FORMULA ###
                    try:
                        formula = self.rpReader.smiles_inchi[meta]['formula']
                    except KeyError:
                        formula = ''
                    ### CHARGE ###
                    try:
                        charge = self.rpReader.model_chemicals[meta]['charge']
                    except KeyError:
                        charge = 0
                    '''
                    #rpsbml.createSpecies(meta, self.chemXref, None, inchi, inchiKey, smiles, compartment_id, charge, formula)
                    #rpsbml.createSpecies(meta, self.chemXref, None, inchi, inchiKey, smiles, compartment_id)
                #4) add the complete reactions and their annotations
                for step in steps:
                    #try:
                        #reac_smiles = self.comp_reac_smiles[step['transformation_id']] ##self.rpReader.rp_transformation[step['transformation_id']]
                    #except KeyError:
                    #    reac_smiles = None
                    #reac_smiles = self.rpReader.rp_transformation[self.rp_paths[pathNum][i+1][comb_path[i]]['transformation_id']]['rule']
                    rpsbml.createReaction('RP'+str(step['step']), # parameter 'name' of the reaction deleted : 'RetroPath_Reaction_'+str(step['step']),
                            'B_999999', #only for genericModel
                            'B_0', #only for genericModel
                            step,
                            self.rpReader.rp_transformation[step['transformation_id']]['rule'],
                            self.reacXref,
                            self.rpReader.rp_transformation,
                            compartment_id)
                #5) adding the consumption of the target
                targetRule = {'rule_id': None, 'left': {[i for i in all_meta if i[:6]=='TARGET'][0]: 1}, 'right': [], 'step': None, 'path_id': None, 'transformation_id': None, 'rule_score': None}
                rpsbml.createReaction('targetSink',
                        'B_999999',
                        'B_0',
                        targetRule,
                        None,
                        self.reacXref,
                        self.rpReader.rp_transformation,
                        compartment_id,
                        True)
                #6) Optional?? Add the flux objectives. Could be in another place, TBD
                #rpsbml.createFluxObj('rpFBA_obj', 'RP0', 1, True)
                rpsbml.createFluxObj('rpFBA_obj', 'targetSink', 1, True)
                #self.sbml_paths['rp_'+str(path_id)] = rpsbml
                self.sbml_paths['rp_'+str(path_id)+'_'+str(altPathNum)] = rpsbml
                altPathNum += 1




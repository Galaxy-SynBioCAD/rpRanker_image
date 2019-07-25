import pickle
import os
import sys
import gzip
import copy
import logging
import itertools
from .rpSBML import rpSBML


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
    def __init__(self):
        ##### stuff to load from cache #####
        self.full_reactions = None
        self.rr_reactions = None
        self.mnxm_strc = None
        self.chemXref = None
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
            self.mnxm_strc = pickle.load(gzip.open(dirname+'/cache/mnxm_strc.pickle.gz', 'rb'))
        except FileNotFoundError as e:
            logging.error(e)
            return False
        try:
            self.chemXref = pickle.load(gzip.open(dirname+'/cache/chemXref.pickle.gz', 'rb'))
        except FileNotFoundError as e:
            logging.error(e)
            return False
        return True

    
    ################################################################
    ######################### PUBLIC FUNCTIONS #####################
    ################################################################


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
            #identify the main compounds and remove them from the full reaction
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
                if i in self.mnxm_strc:
                    if 'smiles' in self.mnxm_strc[i] and not self.mnxm_strc[i]['smiles']==None:
                        reac_smiles = step['reaction_rule'].split('>>')
                        if reac_side=='left':
                            reac_smiles[1] += '.'+self.mnxm_strc[i]['smiles']
                        else: #if it is anything else than left here, above should detect error
                            reac_smiles[0] += '.'+self.mnxm_strc[i]['smiles']
                        step['reaction_rule'] = reac_smiles[0]+'>>'+reac_smiles[1]
                    else:
                        logging.warning('There are no SMILES defined for '+str(i)+' in self.mnxm_strc[i]')
                        continue
                else:
                    logging.warning('Cannot find '+str(i)+' in self.mnxm_strc')
                    continue


    ## Add the cofactors to monocomponent reactions
    #
    # @param step Step in a pathway
    # @param pathway_cmp_mnxm Dictionnary of intermediate compounds with their public ID's
    # @return Boolean determine if the step is to be added
    def addCofactors_step(self, step, pathway_cmp_mnxm):
        try:
            if self.rr_reactions[step['rule_id']]['rel_direction']==-1:
                self.completeReac(step, 'right', self.rr_reactions[step['rule_id']]['left'], self.full_reactions[self.rr_reactions[step['rule_id']]['reac_id']]['right'], pathway_cmp_mnxm)
                self.completeReac(step, 'left', self.rr_reactions[step['rule_id']]['right'], self.full_reactions[self.rr_reactions[step['rule_id']]['reac_id']]['left'], pathway_cmp_mnxm)
            elif self.rr_reactions[step['rule_id']]['rel_direction']==1:
                self.completeReac(step, 'right', self.rr_reactions[step['rule_id']]['left'], self.full_reactions[self.rr_reactions[step['rule_id']]['reac_id']]['left'], pathway_cmp_mnxm)
                self.completeReac(step, 'left', self.rr_reactions[step['rule_id']]['right'], self.full_reactions[self.rr_reactions[step['rule_id']]['reac_id']]['right'], pathway_cmp_mnxm)
            else:
                logging.error('Relative direction can only be 1 or -1: '+str(self.rr_reactions[step['rule_id']]['rel_direction']))
        except KeyError:
            return False
        return True


    ## Function to reconstruct the heterologous pathway
    #
    #  Read each pathway information and RetroRules information to construct heterologous pathways and add the cofactors
    #
    #  @param self Object pointer
    #  @param rpsbml rpSBML object with a single model
    #  @return Boolean if True then you keep that model for the next step, if not then ignore it
    def addCofactors(self, rpsbml, compartment_id='MNXC3'):
        #This keeps the IDs conversions to the pathway
        pathway_cmp_mnxm = {}
        rp_path = rpsbml.outPathsDict()
        ori_rp_path = copy.deepcopy(rp_path)
        #We reverse the loop to ID the intermediate CMP to their original ones
        for stepNum in sorted(list(rp_path), reverse=True):
            if self.addCofactors_step(rp_path[stepNum], pathway_cmp_mnxm):
                ###add the new cofactors to the SBML
                #remove the original species from the monocomponent reaction
                reactants = set(set(rp_path[stepNum]['left'].keys())-set(ori_rp_path[stepNum]['left'].keys()))
                products = set(set(rp_path[stepNum]['right'].keys())-set(ori_rp_path[stepNum]['right'].keys()))
                for species in reactants|products:
                    #check to make sure that they do not yet exist and if not create a new one
                    if not rpsbml.speciesExists(species):
                        xref = {}
                        inchi = None
                        inchikey = None
                        smiles = None
                        try:
                            xref = self.chemXref[species]
                        except KeyError:
                            #TODO: although there should not be any 
                            #intermediate species here consider
                            #removing this warning
                            logging.warning('Cannot find the xref for this species: '+str(species))
                            pass
                        try:
                            inchi = self.mnxm_strc[species]['inchi']
                        except KeyError:
                            logging.warning('Cannot find the inchi for this species: '+str(species))
                            pass
                        try:
                            inchikey = self.mnxm_strc[species]['inchikey']
                        except KeyError:
                            logging.warning('Cannot find the inchikey for this species: '+str(species))
                            pass
                        try:
                            smiles = self.mnxm_strc[species]['smiles']
                        except KeyError:
                            logging.warning('Cannot find the smiles for this species: '+str(species))
                            pass
                        #add the new species to rpsbml
                        rpsbml.createSpecies(species, 
                                xref, 
                                None, 
                                inchi,
                                inchikey,
                                smiles,
                                compartment_id)
                #add the new species to the RP reactions
                reac = rpsbml.model.getReaction(rp_path[stepNum]['reaction_id'])
                for pro in products:
                    prod = reac.createProduct()
                    prod.setSpecies(str(pro)+'__64__'+str(compartment_id))
                    prod.setConstant(True)
                    prod.setStoichiometry(rp_path[stepNum]['right'][pro])
                return True
            else:
                return False

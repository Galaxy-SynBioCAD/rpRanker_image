import pickle
import os
import sys
import gzip
import copy
import itertools
import difflib
from .rpSBML import rpSBML

import logging

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
        self.logger = logging.getLogger(__name__)
        self.logger.debug('Started instance of rpCofactors')
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
            self.deprecatedMNXM_mnxm = pickle.load(open(dirname+'/cache/deprecatedMNXM_mnxm.pickle', 'rb'))
        except FileNotFoundError as e:
            self.logger.error(e)
            return False
        try:
            self.full_reactions = pickle.load(open(dirname+'/cache/full_reactions.pickle', 'rb'))
        except FileNotFoundError as e:
            self.logger.error(e)
            return False
        try:
            self.rr_reactions = pickle.load(open(dirname+'/cache/rr_reactions.pickle', 'rb'))
        except FileNotFoundError as e:
            self.logger.error(e)
            return False
        try:
            self.mnxm_strc = pickle.load(gzip.open(dirname+'/cache/mnxm_strc.pickle.gz', 'rb'))
        except FileNotFoundError as e:
            self.logger.error(e)
            return False
        try:
            self.chemXref = pickle.load(gzip.open(dirname+'/cache/chemXref.pickle.gz', 'rb'))
        except FileNotFoundError as e:
            self.logger.error(e)
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
    def completeReac(self, step, reac_side, f_reac_side, pathway_cmp_mnxm):
        '''
        ## BUG fix, remove MNXM1 (i.e. hydrogen ion) from the rr_reactions since they are commonly not
        # found in the full reaction and causes an error
        hydro_diff = None
        try:
            hydro_diff = {'MNXM1': step[reac_side]['MNXM1']}
            del step[reac_side]['MNXM1']
        except KeyError:
            pass
        '''
        f_reac = self.full_reactions[self.rr_reactions[step['rule_id']][step['rule_mnxr']]['reac_id']][f_reac_side]
        '''
        print('f_reac_side: '+str(f_reac_side))
        print('reac_side: '+str(reac_side))
        print('step[reac_side]: '+str(step[reac_side]))
        print('f_reac: '+str(f_reac))
        '''
        ######## COFACTORS #####
        if reac_side=='right':
            #from the monocomponent side, remove the main species from RR
            noMain_fullReac = {i:f_reac[i] for i in f_reac if i!=self.rr_reactions[step['rule_id']][step['rule_mnxr']]['subs_id']}
            #find the intermediate compound
            cmp_diff = list(step[reac_side].keys()-noMain_fullReac.keys())
            #add it to the dictionnary
            if len(cmp_diff)==1:
                pathway_cmp_mnxm.update({cmp_diff[0]: self.rr_reactions[step['rule_id']][step['rule_mnxr']]['subs_id']})
            else:
                self.logger.warning('There are more than 1 or None of cmp_diff: '+str(cmp_diff))
        elif reac_side=='left':
            #identify the main compounds and remove them from the full reaction
            try:
                #works with out_paths.csv (RP2) but not with RP3... not sure why
                toRem = [pathway_cmp_mnxm[i] for i in list(step[reac_side].keys()-f_reac.keys())]
                noMain_fullReac = {i:f_reac[i] for i in f_reac if i not in toRem}
            except KeyError:
                #INFO1: this happens if the RR reaction does not match the full reaction (original MNX)
                # and usually happens since there can be more than one reaction rule associated
                # with a reaction when the original reaction species do not match
                #INFO2: can also happen when a reaction has more than one product that need to be recovered here
                # Means that there are two species that need to be recovered instead of the one from
                # excpected monocomponent reaction
                #Temp fix: use the full_reactions main species to do it
                #works with RP3, but don't like to need to refer to the original reaction
                try:
                    self.logger.debug('Reverting to using the full reactions main reaction')
                    ## control for INFO2
                    toRem = [i for i in list(step[reac_side].keys()-f_reac.keys())]
                    if len(toRem)==1:
                        main = self.full_reactions[self.rr_reactions[step['rule_id']][step['rule_mnxr']]['reac_id']]['main_'+str(f_reac_side)]
                        noMain_fullReac = {i:f_reac[i] for i in f_reac if i not in main}
                        pathway_cmp_mnxm.update({list(step[reac_side].keys()-f_reac.keys())[0]: main[0]})
                    else:
                        self.logger.warning('There are more than one unknown species. Attempting a rescue...')
                        main = self.full_reactions[self.rr_reactions[step['rule_id']][step['rule_mnxr']]['reac_id']]['main_'+str(f_reac_side)]
                        unknownMain = list(step[reac_side].keys()-f_reac.keys())[0]
                        pathway_cmp_mnxm.update({unknownMain: main[0]})
                        #recover the inchi keys from species in the full reactions, remove the knowns,
                        toCompare = [i for i in f_reac if not i==unknownMain]
                        inchikey_freac = {}
                        for mnxm in f_reac:
                            try:
                                inchikey_freac[self.mnxm_strc[mnxm]['inchikey'].replace('-','')] = mnxm
                            except KeyError:
                                self.logger.error('Cannot find the inchikey for '+str(mnxm))
                                return False
                        toFind = [i for i in toRem if not i==main[0]]
                        for i in toFind:
                            inchikey = difflib.get_close_matches(i, list(inchikey_freac.keys()), n=1)[0]
                            #main.append(inchikey_freac[inchikey])
                            #pathway_cmp_mnxm.update({inchikey: inchikey_freac[inchikey]})
                        noMain_fullReac = {i:f_reac[i] for i in f_reac if i not in main}
                    '''
                    print('#######################')
                    print(toRem)
                    print(main)
                    print(noMain_fullReac)
                    print(step[reac_side].keys())
                    print('#######################')
                    '''
                except (KeyError, IndexError) as e:
                    self.logger.error('Could not find intermediate in pathway_cmp_mnxm '+str(e))
                    return False
        else:
            self.logger.warning('Direction can only be right or left')
            return False
        #calculate the difference between the two
        diff = {i: noMain_fullReac[i] for i in noMain_fullReac.keys()-step[reac_side].keys()}
        #update the reaction
        step[reac_side].update(diff)
        '''
        if not hydro_diff==None:
            step[reac_side].update(hydro_diff)
        '''
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
                    self.logger.warning('Could not find the intermediate compound in full reaction: '+str(i))
                    self.logger.warning('Setting to 1.0')
                    step[reac_side][i] = 1.0
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
                        #TODO: if any of the steps fail you should revert to the original reaction rule
                        self.logger.warning('There are no SMILES defined for '+str(i)+' in self.mnxm_strc')
                        continue
                else:
                    self.logger.debug('Cannot find '+str(i)+' in self.mnxm_strc')
                    continue
        return True

    ## Add the cofactors to monocomponent reactions
    #
    # @param step Step in a pathway
    # @param pathway_cmp_mnxm Dictionnary of intermediate compounds with their public ID's
    # @return Boolean determine if the step is to be added
    def addCofactors_step(self, step, pathway_cmp_mnxm):
        '''
        print('########## '+str(step['rule_mnxr'])+' ############')
        print(pathway_cmp_mnxm)
        print(step['left'])
        print(step['right'])
        print('----------------------')
        '''
        if self.rr_reactions[step['rule_id']][step['rule_mnxr']]['rel_direction']==-1:
            if not self.completeReac(step, 'right', 'right', pathway_cmp_mnxm):
                self.logger.error('Could not recognise reaction rule for step {}'.format(step))
                return False
            if not self.completeReac(step, 'left', 'left', pathway_cmp_mnxm):
                self.logger.error('Could not recognise reaction rule for step {}'.format(step))
                return False
        elif self.rr_reactions[step['rule_id']][step['rule_mnxr']]['rel_direction']==1:
            if not self.completeReac(step, 'right', 'left', pathway_cmp_mnxm):
                self.logger.error('Could not recognise reaction rule for step {}'.format(step))
                return False
            if not self.completeReac(step, 'left', 'right', pathway_cmp_mnxm):
                self.logger.error('Could not recognise reaction rule for step {}'.format(step))
                return False
        else:
            self.logger.error('Relative direction can only be 1 or -1: '+str(self.rr_reactions[step['rule_id']][step['rule_mnxr']]['rel_direction']))
            return False
        '''
        print('----------------------')
        print(step['left'])
        print(step['right'])
        print(pathway_cmp_mnxm)
        print('----------------------')
        '''
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
                            try:
                                xref = self.chemXref[self.deprecatedMNXM_mnxm[species]]
                            except KeyError:
                                #TODO: although there should not be any
                                #intermediate species here consider
                                #removing this warning
                                self.logger.warning('Cannot find the xref for this species: '+str(species))
                                pass
                        try:
                            inchi = self.mnxm_strc[species]['inchi']
                        except KeyError:
                            try:
                                inchi = self.mnxm_strc[self.deprecatedMNXM_mnxm[species]]['inchi']
                            except KeyError:
                                self.logger.warning('Cannot find the inchi for this species: '+str(species))
                                pass
                        try:
                            inchikey = self.mnxm_strc[species]['inchikey']
                        except KeyError:
                            try:
                                inchikey = self.mnxm_strc[self.deprecatedMNXM_mnxm[species]]['inchikey']
                            except KeyError:
                                self.logger.warning('Cannot find the inchikey for this species: '+str(species))
                                pass
                        try:
                            smiles = self.mnxm_strc[species]['smiles']
                        except KeyError:
                            try:
                                smiles = self.mnxm_strc[self.deprecatedMNXM_mnxm[species]]['smiles']
                            except KeyError:
                                self.logger.warning('Cannot find the smiles for this species: '+str(species))
                                pass
                        #add the new species to rpsbml
                        rpsbml.createSpecies(species,
                                compartment_id,
                                xref,
                                None,
                                inchi,
                                inchikey,
                                smiles)
                #add the new species to the RP reactions
                reac = rpsbml.model.getReaction(rp_path[stepNum]['reaction_id'])
                for pro in products:
                    prod = reac.createProduct()
                    prod.setSpecies(str(pro)+'__64__'+str(compartment_id))
                    prod.setConstant(True)
                    prod.setStoichiometry(rp_path[stepNum]['right'][pro])
                for sub in reactants:
                    subs = reac.createReactant()
                    subs.setSpecies(str(sub)+'__64__'+str(compartment_id))
                    subs.setConstant(True)
                    subs.setStoichiometry(rp_path[stepNum]['left'][sub])
            else:
                #if the cofactors cannot be found delete it from the list
                return False
        return True

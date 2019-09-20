import pickle
import os
import sys
import gzip
import copy
import itertools
#import difflib
import rpSBML

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
        self.logger.info('Started instance of rpCofactors')
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


    ## Private function to fetch the required data, parse them and generate the pickle
    #
    #  Opens the previously generated cache to the object memory
    #
    # @param The oject pointer
    # @return Boolean detemining the success of the function or not
    def _loadCache(self, fetchInputFiles=False):
        dirname = os.path.dirname(os.path.abspath( __file__ ))
        #################### make the local folders ############################
        # input_cache
        if not os.path.isdir(dirname+'/input_cache')
            os.mkdir(dirname+'/input_cache')
        # cache
        if not os.path.isdir(dirname+'/cache'):
            os.mkdir(dirname+'/cache')
        ###################### Fetch the files if necessary ######################
        #chem_xref
        if not os.path.isfile(dirname+'/input_cache/chem_xref.tsv') or fetchInputFiles:
            urllib.request.urlretrieve('https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv', 
                    dirname+'/input_cache/chem_xref.tsv')
        # rr_compounds.tsv
        #TODO: need to add this file to the git or another location
        if not os.path.isfile(dirname+'/input_cache/rr_compounds.tsv') or fetchInputFiles:
            urllib.request.urlretrieve(
                    'TOADD', 
                    dirname+'/input_cache')
            #tf = tarfile.open(dirname+'/input_cache/retrorules_preparsed.tar.xz')
            #tf.extractall(path=dirname+'/input_cache/')
            #tf.close()
        # rules_rall.tsv
        if not os.path.isfile(dirname+'/input_cache/rules_rall.tsv') or fetchInputFiles:
            urllib.request.urlretrieve('TOADD', 
                    dirname+'/input_cache/')
        # rxn_recipes.tsv
        if not os.path.isfile(dirname+'/input_cache/rxn_recipes.tsv') or fetchInputFiles:
            urllib.request.urlretrieve('TOADD', 
                    dirname+'/input_cache/')
        # chem_prop.tsv
        if not os.path.isfile(dirname+'/input_cache/chem_prop.tsv') or fetchInputFiles:
            urllib.request.urlretrieve('https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv', 
                    dirname+'/input_cache/chem_prop.tsv')
        ###################### Populate the cache #################################
        if not os.path.isfile(dirname+'/cache/deprecatedMNXM_mnxm.pickle'):
            rpCache.deprecatedMNXM(dirname+'/input_cache/chem_xref.tsv')
            pickle.dump(rpCache.deprecatedMNXM_mnxm, open(dirname+'/cache/deprecatedMNXM_mnxm.pickle', 'wb'))
        self.deprecatedMNXM_mnxm = pickle.load(open(dirname+'/cache/deprecatedMNXM_mnxm.pickle', 'rb'))
        if not os.path.isfile(dirname+'/cache/mnxm_strc.pickle.gz'):
            pickle.dump(rpCache.mnx_strc(dirname+'/input_cache/rr_compounds.tsv', 
                            dirname+'/input_cache/chem_prop.tsv'), 
                        gzip.open(dirname+'/cache/mnxm_strc.pickle.gz','wb'))
        self.mnxm_strc = pickle.load(gzip.open(dirname+'/cache/mnxm_strc.pickle.gz', 'rb'))
        if not os.path.isfile(dirname+'/cache/full_reactions.pickle'):
            pickle.dump(rpCache.full_reac('input_cache/rxn_recipes.tsv'),
                    open('cache/full_reactions.pickle', 'wb'))
        self.full_reactions = pickle.load(open(dirname+'/cache/full_reactions.pickle', 'rb'))
        if not os.path.isfile(dirname+'/cache/chemXref.pickle.gz'):
            pickle.dump(rpCache.mnx_chemXref(dirname+'/input_cache/chem_xref.tsv'), 
                    gzip.open(dirname+'/cache/chemXref.pickle.gz','wb'))
        self.chemXref = pickle.load(gzip.open(dirname+'/cache/chemXref.pickle.gz', 'rb'))
        if not os.path.isfile(dirname+'/cache/rr_reactions.pickle'):
            pickle.dump(
                    rpCache.retro_reactions(dirname+'/input_cache/rules_rall.tsv'), 
                    open(dirname+'/cache/rr_reactions.pickle', 'wb'))
        self.rr_reactions = pickle.load(open(dirname+'/cache/rr_reactions.pickle', 'rb'))
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
    def completeReac(self, step, rr_reac, full_reac, mono_side, rr_string, pathway_cmp_mnxm):
        if mono_side:
            ## add the unknown species to pathway_cmp_mnxm for the next steps
            rr_mono_cmp = list(rr_reac.keys())
            step_mono_cmp = list(step.keys())
            if (len(rr_mono_cmp)==1 and len(step_mono_cmp)==1):
                #this is purposely overwitten since the main cmp between reactions can change
                pathway_cmp_mnxm[step_mono_cmp[0]] = rr_mono_cmp[0]
            else:
                self.logger.warning('There should be only one compound on the left for monocomponent reaction: rr_mono_cmp: '+str(rr_mono_cmp)+' step_mono_cmp: '+str(step_mono_cmp))
                return False
        ## add the side species
        for toAdd in full_reac.keys()-rr_reac.keys():
            step.update({toAdd: full_reac[toAdd]})
            ### update the reaction rule string
            try:
                smi = self.mnxm_strc[toAdd]['smiles']
                if not smi==None:
                    rr_string += '.'+str(smi)
            except KeyError:
                self.logger.warning('Cannot find smiles structure for '+str(toAdd))
        ## Update the the stochio
        for step_spe in step:
            if step_spe in full_reac:
                if not step[step_spe]==full_reac[step_spe]:
                    stochio_diff = full_reac[step_spe]-step[step_spe]
                    step[step_spe] = full_reac[step_spe]
                    if stochio_diff<0:
                        self.logger.warning('full_reac stochio should never be smaller than step')
                        continue
                    for i in range(stochio_diff):
                        ### update the reaction rule string
                        try:
                            smi = self.mnxm_strc[step_spe]['smiles']
                            if not smi==None:
                                rr_string += '.'+str(smi)
                        except KeyError:
                            self.logger.warning('Cannot find smiles structure for '+str(toAdd))
            elif step_spe in pathway_cmp_mnxm:
                if pathway_cmp_mnxm[step_spe] in full_reac:
                    if not step[step_spe]==full_reac[pathway_cmp_mnxm[step_spe]]:
                        step[step_spe] = full_reac[pathway_cmp_mnxm[step_spe]]
            #Its fine if the stochio is not updated, better than ignoring a whole pathway
                #else:
                #    self.logger.warning('Cannot find '+str(step_spe)+' in full reaction')
                #    return False
            #else:
            #    self.logger.warning('Cannot find '+str(step_spe)+' in pathway_cmp_mnxm')
            #    return False
        return True, rr_string


    ## Add the cofactors to monocomponent reactions
    #
    # @param step Step in a pathway
    # @param pathway_cmp_mnxm Dictionnary of intermediate compounds with their public ID's
    # @return Boolean determine if the step is to be added
    def addCofactors_step(self, step, pathway_cmp_mnxm):
        reac_smiles_left = step['reaction_rule'].split('>>')[0]
        reac_smiles_right = step['reaction_rule'].split('>>')[1]
        if self.rr_reactions[step['rule_id']][step['rule_mnxr']]['rel_direction']==-1:
            isSuccess, reac_smiles_right = self.completeReac(step['right'], 
                    self.rr_reactions[step['rule_id']][step['rule_mnxr']]['left'],
                    self.full_reactions[step['rule_mnxr']]['right'], 
                    True,
                    reac_smiles_right,
                    pathway_cmp_mnxm)
            if not isSuccess:
                self.logger.error('Could not recognise reaction rule for step '+str(step))
                return False
            isSuccess, reac_smiles_left = self.completeReac(step['left'],
                    self.rr_reactions[step['rule_id']][step['rule_mnxr']]['right'],
                    self.full_reactions[step['rule_mnxr']]['left'], 
                    False,
                    reac_smiles_left,
                    pathway_cmp_mnxm)
            if not isSuccess:
                self.logger.error('Could not recognise reaction rule for step '+str(step))
                return False
        elif self.rr_reactions[step['rule_id']][step['rule_mnxr']]['rel_direction']==1:
            isSuccess, reac_smiles_right = self.completeReac(step['right'],
                    self.rr_reactions[step['rule_id']][step['rule_mnxr']]['left'],
                    self.full_reactions[step['rule_mnxr']]['left'], 
                    True,
                    reac_smiles_right,
                    pathway_cmp_mnxm)
            if not isSuccess:
                self.logger.error('Could not recognise reaction rule for step '+str(step))
                return False
            isSuccess, reac_smiles_left = self.completeReac(step['left'],
                    self.rr_reactions[step['rule_id']][step['rule_mnxr']]['right'],
                    self.full_reactions[step['rule_mnxr']]['right'], 
                    False,
                    reac_smiles_left,
                    pathway_cmp_mnxm)
            if not isSuccess:
                self.logger.error('Could not recognise reaction rule for step '+str(step))
                return False
        else:
            self.logger.error('Relative direction can only be 1 or -1: '+str(self.rr_reactions[step['rule_id']][step['rule_mnxr']]['rel_direction']))
            return False
        step['reaction_rule'] = reac_smiles_left+'>>'+reac_smiles_right
        return True


    ## Function to reconstruct the heterologous pathway
    #
    #  Read each pathway information and RetroRules information to construct heterologous pathways and add the cofactors
    #
    #  @param self Object pointer
    #  @param rpsbml rpSBML object with a single model
    #  @return Boolean if True then you keep that model for the next step, if not then ignore it
    def addCofactors(self, rpsbml, compartment_id='MNXC3', pathId='rp_pathway'):
        #This keeps the IDs conversions to the pathway
        pathway_cmp_mnxm = {}
        rp_path = rpsbml.outPathsDict(pathId)
        ori_rp_path = copy.deepcopy(rp_path)
        #We reverse the loop to ID the intermediate CMP to their original ones
        for stepNum in sorted(list(rp_path), reverse=True):
        #for stepNum in sorted(list(rp_path)):
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
                        try:
                            chemName = self.mnxm_strc[species]
                        except KeyError:
                            chemName = None
                        rpsbml.createSpecies(species,
                                compartment_id,
                                chemName,
                                xref,
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

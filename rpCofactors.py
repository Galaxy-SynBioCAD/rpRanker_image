import pickle
import os
import gzip
import copy
from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi, MolToInchiKey, AddHs
import logging

import rpSBML
import rpThermo

import sys
import gc

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
    def __init__(self, rpReader, userXrefDbName=None):
        self.comp_reac_smiles = {}
        self.rpReader = rpReader
        self.chem_xref = None
        self.convertid_inchi = {}
        self.sbml_paths = {}
        self.full_reactions = None
        self.userXrefDbName = userXrefDbName #TODO: change to interpret all inputs xref database from the input file directly 
        self.rr_reactions = None
        if not self._loadCache(os.getcwd()+'/cache'):
            raise ValueError

    ##
    #
    #
    def _checkFilePath(self, path, filename):
        """Check that the directory and the filename are valid and choose to use
        either the local or the global path
        """
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


    ## Load the cache required for this class
    #
    #
    def _loadCache(self, path):
        try:
        	self.chem_xref = pickle.load(gzip.open(self._checkFilePath(path, 'chemXref.pickle.gz'), 'rb'))
        except FileNotFoundError:
            logging.error('The file '+str(path+'/chemXref.pickle.gz')+' does not seem to exists')
            return False
        try:
            self.full_reactions = pickle.load(open(self._checkFilePath(path, 'full_reactions.pickle'), 'rb'))
        except FileNotFoundError:
            logging.error('The file '+str(path, '/full_reactions.pickle')+' does not seem to exist')
            return False
        try:
            self.rr_reactions = pickle.load(open(self._checkFilePath(path, 'rr_reactions.pickle'), 'rb'))
        except FileNotFoundError:
            logging.error('The file '+str(path, '/rr_reactions.pickle')+' does not seem to exist')
            return False
        return True

    ## Function to identify the user identifier from the MNX one
    #
    #  Go through the cross references files to extract the corresponding identifiers
    #
    #  @param self Object pointer
    #  @param compound The MNX identifier to be identifiy
    #  @param db The user database containing a cross reference
    #  @return finalID 
    def cmpd_identification_xref(self, compound, db):
        if not 'MNXM' in compound[:4]:
            if compound in self.rpReader.in_xref:
                db_cid = self.rpReader.in_xref[compound][db]
                if db_cid in self.chem_xref[db]:  #for i in self.chem_xref if db in self.chem_xref[i]: for n in self.chem_xref[i][db]: if db_cid == n:
                    final_id = self.chem_xref[db][db_cid]
        else:
            tmp_id = self.chem_xref[compound][db]
            for i in self.in_xref:
                for n in tmp_id:
                    if n == self.rpReader.in_xref[i][db]:
                        final_id = i
        return final_id


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
    def convert_depiction(self, idepic, itype='smiles', otype={'inchikey'}):
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
            for i in self.rpReader.in_inchi:
                c = self.convert_depiction(idepic=self.rpReader.in_inchi[i], itype='inchi', otype={'inchikey'})
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


    ## Function to reconstruct the heterologous pathway
    #
    #  Read each pathway information and RetroRules information to construct heterologous pathways and add the cofactors
    #
    #  @param self Object pointer
    #  @param in_rp_paths All pathways
    #  @param rr_reactions Dictionnary of reactions ID (rules_rall)
    #  @param full_reactions Dictionnary of reactions of origin (rxn_recipes)
    #  @param rp_smiles Dictionnary of smile and structure for each compound (compounds.txt)
    #  @param self.rpReader.rp_transformation Dictionnary discribing each transformation
    #  @param smiles_inchi dictionnary describing the inchi for each smile (scope)
    #  @param rp_smiles_inchi dictionnary describing the inchi for each smile
    #  @return out_rp_paths Complete heterologous pathway object
    def addCofactors(self):
        self.comp_reac_smiles = {}
        rp_paths = copy.deepcopy(self.rpReader.rp_paths)
        for path in rp_paths:
            for step in path:
                try:
                    reac = copy.deepcopy(self.full_reactions[self.rr_reactions[step['rule_id']]['reaction']] )
                except KeyError:
                    logging.error('Cannot find rule_id: '+step['rule_id'])
                    sys.exit('Cannot find rule_id. You are probably not using the last configuration of retro_rules')
                    break
                rp_path_left = None
                rp_path_right = None
                #use the relative direction in retro_rules to determine the direction of the origin reaction
                if self.rr_reactions[step['rule_id']]['rel_direction'] == '-1':
                    rp_path_right = 'right'
                    rp_path_left = 'left'
                if self.rr_reactions[step['rule_id']]['rel_direction'] == '1':
                    rp_path_right = 'left'
                    rp_path_left = 'right'
                if rp_path_left==None and rp_path_right==None:
                    logging.error('Cannot find the direction from the rule compared to rule_id')
                    break
                ##remove the main right and main left compounds from the full reaction
                ##############RIGHT######################
                toremove_l = []
                toremove_r = []
                for i in (self.rr_reactions[step['rule_id']][rp_path_right]).split('.'):
                    if self.rr_reactions[step['rule_id']]['rel_direction'] == '-1':
                        for n in reac[rp_path_left].keys():
                            if i == n and not i in toremove_l:
                                toremove_l.append(i)
                    elif self.rr_reactions[step['rule_id']]['rel_direction'] == '1':
                        for n in reac[rp_path_right].keys():
                            if i == n and not i in toremove_r:
                                toremove_r.append(i)
                for i in toremove_l:
                    del reac[rp_path_left][i]
                for i in toremove_r:
                    del reac[rp_path_right][i]
                ##############LEFT########################
                toremove_l = []
                toremove_r = []
                for i in (self.rr_reactions[step['rule_id']][rp_path_left]).split('.'):
                    if self.rr_reactions[step['rule_id']]['rel_direction'] == '-1':
                        for n in reac[rp_path_right].keys():
                            if i == n and not i in toremove_r:
                                toremove_r.append(i)
                    elif self.rr_reactions[step['rule_id']]['rel_direction'] == '1':
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
                    while toAdd in self.rpReader.deprecatedMNXM_mnxm:
                        toAdd = self.rpReader.deprecatedMNXM_mnxm[toAdd]
                    if not toAdd in step['left']:
                        if ('MNXM' in toAdd[:4]) and (self.rpReader.in_xref is not None or self.rpReader.in_inchi is not None):
                            try:
                                tmp_id = self.cmpd_identification_xref(toAdd, self.userXrefDbName)
                                new_toAdd = tmp_id
                                step['left'][new_toAdd+':'+toAdd] = toAdd_left[toAdd]
                            except(KeyError, UnboundLocalError, TypeError):
                                try:
                                    tmp_id = self.cmpd_identification_inchi(toAdd)
                                    new_toAdd = tmp_id
                                    step['left'][new_toAdd+':'+toAdd] = toAdd_left[toAdd]
                                except(KeyError, UnboundLocalError):
                                    logging.warning("No identifier have been found for the compound "+str(toAdd))
                                    step['left'][toAdd] = toAdd_left[toAdd]
                                except TypeError:
                                    logging.warning("You did not provide a sink file")
                for toAdd in toAdd_right:
                    while toAdd in self.rpReader.deprecatedMNXM_mnxm:
                        toAdd = self.rpReader.deprecatedMNXM_mnxm[toAdd]
                    if not toAdd in step['right']:
                        if ('MNXM' in toAdd[:4]) and (self.rpReader.in_xref is not None or self.rpReader.in_inchi is not None):
                            try:
                                tmp_id = self.cmpd_identification_xref(toAdd, self.userXrefDbName)
                                new_toAdd = tmp_id
                                step['right'][new_toAdd+':'+toAdd] = toAdd_right[toAdd]
                            except(KeyError, UnboundLocalError, TypeError):
                                try:
                                    tmp_id = self.cmpd_identification_inchi(toAdd)
                                    new_toAdd = tmp_id
                                    step['right'][new_toAdd+':'+toAdd] = toAdd_right[toAdd]
                                except(KeyError, UnboundLocalError):
                                    logging.warning("No identifier have been found for the compound "+str(toAdd))
                                    step['right'][toAdd] = toAdd_right[toAdd]
                                except TypeError:
                                    logging.warning("You did not provide a sink file")
                #reconstruct the complete reaction smiles by adding the smiles of cofactors
                if not step['transformation_id'] in self.comp_reac_smiles:
                    try:
                        scope_smiles = self.rpReader.rp_transformation[step['transformation_id']].split('>>')
                        tmp_smiles_l=''
                        tmp_smiles_r=''
                        if self.rr_reactions[step['rule_id']]['rel_direction'] == '-1':
                            smiles_l = scope_smiles[0]
                            smiles_r = scope_smiles[1]
                            adding_left = toAdd_right
                            adding_right = toAdd_left
                        if self.rr_reactions[step['rule_id']]['rel_direction'] == '1':
                            smiles_l = scope_smiles[1]
                            smiles_r = scope_smiles[0]
                            adding_left = toAdd_left
                            adding_right = toAdd_right
                        for toAdd in adding_left:
                            tmp = {}
                            try:
                                #TODO: Hack, need to get to the source of this
                                if self.rpReader.rp_smiles_inchi[self.rpReader.rp_smiles[toAdd]['smiles']]=='InChI=1S':
                                    tmp['smiles'] = '[H]' 
                                else:
                                    tmp = self.convert_depiction(idepic=self.rpReader.rp_smiles_inchi[self.rpReader.rp_smiles[toAdd]['smiles']], itype='inchi', otype={'smiles'})
                                i = 0
                                while i < int(adding_left[toAdd]):
                                    tmp_smiles_l += tmp['smiles']+'.'
                                    i += 1
                            except KeyError:
                                try:
                                    #TODO: hack, need to get to the source of this
                                    if self.rpReader.smiles_inchi[toAdd]['inchi'][0:5]=='InChI':
                                        tmp = self.convert_depiction(idepic=self.rpReader.smiles_inchi[toAdd]['inchi'], itype='inchi', otype={'smiles'})
                                    i = 0
                                    while i < int(adding_left[toAdd]):
                                        tmp_smiles_l += tmp['smiles']+'.'
                                        i += 1
                                except KeyError:
                                    tmp_smiles_l=''
                        tmp_smiles_l += smiles_l
                        for toAdd in adding_right:
                            tmp = {}
                            try:
                                #TODO: Hack, need to ge to the source of this
                                if self.rpReader.rp_smiles_inchi[self.rpReader.rp_smiles[toAdd]['smiles']]=='InChI=1S':
                                    tmp['smiles'] = '[H]' 
                                else:
                                    tmp = self.convert_depiction(idepic=self.rpReader.rp_smiles_inchi[self.rpReader.rp_smiles[toAdd]['smiles']], itype='inchi', otype={'smiles'})
                                i = 0
                                while i < int(adding_right[toAdd]):
                                    tmp_smiles_r += tmp['smiles']+'.'
                                    i += 1
                            except KeyError:
                                try:
                                    #TODO: hack, need to get to the source of this
                                    if self.rpReader.smiles_inchi[toAdd]['inchi'][0:5]=='InChI':
                                        tmp = self.convert_depiction(idepic=self.rpReader.smiles_inchi[toAdd]['inchi'], itype='inchi', otype={'smiles'})
                                    i = 0
                                    while i < int(adding_right[toAdd]):
                                        tmp_smiles_r += tmp['smiles']+'.'
                                        i += 1
                                except KeyError:
                                    tmp_smiles_r=''
                        tmp_smiles_r += smiles_r
                        if self.rr_reactions[step['rule_id']]['rel_direction'] == '-1':
                            self.comp_reac_smiles[step['transformation_id']] = tmp_smiles_r+str('>>')+tmp_smiles_l
                        elif self.rr_reactions[step['rule_id']]['rel_direction'] == '1':
                            self.comp_reac_smiles[step['transformation_id']] = tmp_smiles_l+str('>>')+tmp_smiles_r
                    except KeyError:
                        pass
        return rp_paths

    def _get_obj_size(self, obj):
        marked = {id(obj)}
        obj_q = [obj]
        sz = 0

        while obj_q:
            sz += sum(map(sys.getsizeof, obj_q))

            # Lookup all the object reffered to by the object in obj_q.
            # See: https://docs.python.org/3.7/library/gc.html#gc.get_referents
            all_refr = ((id(o), o) for o in gc.get_referents(*obj_q))

            # Filter object that are already marked.
            # Using dict notation will prevent repeated objects.
            new_refr = {o_id: o for o_id, o in all_refr if o_id not in marked and not isinstance(o, type)}

            # The new obj_q will be the ones that were not marked,
            # and we will update marked with their ids so we will
            # not traverse them again.
            obj_q = new_refr.values()
            marked.update(new_refr.keys())

        return sz


    def pathsToRPsbml(self, outputFolderName='sbml_models', compartment_id='MNXC3'):
        #if the output folder does not exist then create it
        #for path in rp_paths:
        sbmlThermo = rpThermo.rpThermo()
        if not os.path.exists(os.getcwd()+'/'+outputFolderName):
            os.makedirs(os.getcwd()+'/'+outputFolderName)
        for path in self.addCofactors():
            steps = [i for i in path]
            path_id = steps[0]['path_id']
            #logging.info('############ Generating an sbml object for path '+str(path_id)+' ###########')
            ###### create a libSBML model ####
            rpsbml = rpSBML.rpSBML()
            #1) create a generic Model, ie the structure and unit definitions that we will use the most
            rpsbml.genericModel('RetroPath_Pathway_'+str(path_id), 'RP_model'+str(path_id))
            #2) create the pathway (groups)
            rpsbml.createPathway('hetero_pathway')
            #3) find all the unique species and add them to the model
            all_meta = set([i for step in steps for lr in ['left', 'right'] for i in step[lr]])
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
                ### INCHI ### 
                try:
                	inchi = self.rpReader.smiles_inchi[meta]['inchi']
                except KeyError:
                    try:
                        inchi = self.rpReader.rp_smiles_inchi[self.rpReader.rp_smiles[meta]['smiles']]
                    except KeyError:
                        inchi = None
                ### SMILES ###
                try:
                    smiles = self.rpReader.rp_smiles[meta]['smiles']
                except KeyError:
                    try:
                        smiles = self.rpReader.smiles_inchi[meta]['smiles']
                    except KeyError:
                        smiles = None
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
                ### COMPOUND IDENTIFICATION ###
                #if 'TARGET' in meta[:6]:  ##'CMPD' in meta[:4] or 'TARGET' in meta[:6]
                #    tmp_meta = cmpd_meta[1]
                #else:
                #    tmp_meta = cmpd_meta[0]
                try:
                    tmp_xref = {}
                    if cmpd_meta[1]:
                        for db in self.chem_xref[cmpd_meta[1]]:
                            tmp_xref[db] = self.chem_xref[cmpd_meta[1]][db]
                    rpsbml.createSpecies(meta, tmp_xref, None, inchi, smiles, compartment_id, charge, formula)
                except (KeyError, IndexError):
                    rpsbml.createSpecies(meta, None, None, inchi, smiles, compartment_id, charge, formula)
            #4) add the complete reactions and their annotations
            for step in path:
                try:
                    reac_smiles = self.comp_reac_smiles[step['transformation_id']] ##self.rpReader.rp_transformation[step['transformation_id']]
                except KeyError:
                    reac_smiles = None
                rpsbml.createReaction('RP'+str(step['step']), # parameter 'name' of the reaction deleted : 'RetroPath_Reaction_'+str(step['step']),
                        'B_INF', #only for genericModel
                        'B_0', #only for genericModel
                        step,
                        reac_smiles,
                        compartment_id)
            '''#5) adding the export reaction
            rpsbml.createCompartment(1, 'e', 'extracellular')
            target_name = [i for i in all_meta if i[:6]=='TARGET'][0]
            rpsbml.createSpecies(target_name, None, None, None, None, 'e')
            for step in path:
                print(step)
                for i in step['right']:
                    if 'TARGET' in i[:6]:
                        rpsbml.createReaction('Target_sink', 'B_0', 'B__999999', step, None, 'e')'''
            #6) Optional?? Add the flux objectives. Could be in another place, TBD
            rpsbml.createFluxObj('rpFBA_obj', 'RP0', 1, True)
            sbmlThermo.pathway_drG_prime_m(rpsbml)
            #print(sys.getsizeof(rpsbml))
            if 100000000<self._get_obj_size(self.sbml_paths):
                print('TOO LLAAAARGE')
                sys.exit()
            print(self._get_obj_size(self.sbml_paths))
            #rpsbml.writeSBML('RP'+str(path_id), os.getcwd()+'/'+outputFolderName)
            self.sbml_paths['rp_'+str(path_id)] = rpsbml

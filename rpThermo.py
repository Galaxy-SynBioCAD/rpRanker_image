import numpy as np
import logging
import pybel
import json
from scipy.special import logsumexp
import copy

#local package
import component_contribution


class Thermodynamics:
    """Combination of equilibrator and group_contribution analysis to calculate the thermodymaics of the individual
    metabolic pathways output from RP2paths and thus RetroPath2.0
    """
    def __init__(self, mnxm_dG, cc_preprocess, deprecatedMNXM_mnxm, pH=7.0, pMg=14.0, I=0.1, temperature=298.15):
        #self.compound_dict = compound_dict
        self.mnxm_dG = mnxm_dG
        self.cc_preprocess = cc_preprocess
        self.deprecatedMNXM_mnxm = deprecatedMNXM_mnxm
        #### set the physiological parameters
        self.pH = pH
        self.pMg = pMg
        self.I = I
        #Constants
        self.R = 8.31e-3   # kJ/(K*mol)
        self.temperature = temperature  # K
        self.phase = 'aqueous'
        self.RT = self.R*self.temperature
        self.RTlog10 = self.RT*np.log(10)
        self.debye_hueckle_a = 2.91482
        self.debye_hueckle_b = 1.6
        self.mg_formation_energy = -455.3  # kJ/mol, formation energy of Mg2+
        self.faraday = 96.485  # kC/mol
        self.conc_lb = 1e-6
        self.conc_ub = 1e-2
        ##### component contribution
        self.min_pH = 0.0
        self.max_pH = 14.0
        #TODO: test
        self.groups_data = component_contribution.inchi2gv.init_groups_data()
        self.group_names = self.groups_data.GetGroupNames()
        self.decomposer = component_contribution.inchi2gv.InChIDecomposer(self.groups_data)
        # Approximation of the temperature dependency of ionic strength effects
        self.DH_alpha = 1e-3*(9.20483*self.temperature) - 1e-5*(1.284668 * self.temperature**2) + 1e-8*(4.95199 * self.temperature**3)
        self.DH_beta = 1.6
        # Debye-Huckel
        self.debye_huckel = self.DH_alpha * self.I **(0.5) / (1.0 + self.DH_beta * self.I**(0.5))
        #temporarely save the calculated dG's
        self.calculated_dG = {}

    ################ PRIVATE FUNCTIONS ###################


    def _select_mnxm_dG(self, mnxm):
        """Given that there can be multiple precalculated dG for a given molecule (MNXM)
        this function gives the priority to a particular order for a given MNXM ID
        alberty>smallest(KEGG_ID)>other(KEGG_ID)
        """
        if mnxm in self.mnxm_dG:
            if 'component_contribution' in self.mnxm_dG[mnxm] and self.mnxm_dG[mnxm]['component_contribution']:
                ####### select the smallest one
                '''
                if len(self.mnxm_dG[mnxm]['component_contribution'])==1:
                    return self.mnxm_dG[mnxm]['component_contribution'][0]['compound_index'], self.mnxm_dG[mnxm]['component_contribution'][0]['group_vector'], self.mnxm_dG[mnxm]['component_contribution'][0]['pmap']['species']
                else:
                    #select the lowest CID and make sure that there
                    toRet = {'CID': 'C99999'}
                    for cmp_dict in self.mnxm_dG[mnxm]['component_contribution']:
                        if int(cmp_dict['CID'][1:])<int(toRet['CID'][1:]):
                            toRet = cmp_dict
                    return toRet['compound_index'], toRet['group_vector'], toRet['pmap']['species']
                '''
                ###### select the one with the most information
                #return the smallest one that has all the information required
                for cmp_d in sorted(self.mnxm_dG[mnxm]['component_contribution'], key=lambda k: int(k['CID'][1:])):
                    if 'compound_index' in cmp_d and 'group_vector' in cmp_d and 'pmap' in cmp_d and 'species' in cmp_d['pmap']:
                        return cmp_d['compound_index'], cmp_d['group_vector'], cmp_d['pmap']['species']
                #if cannot find return the smallest one and return the info that you can
                compound_index = None
                group_vector = None
                pmap_species = None
                try:
                    compound_index = self.mnxm_dG[mnxm]['component_contribution'][0]['compound_index'] 
                except KeyError:
                    pass
                try:
                    group_vector = self.mnxm_dG[mnxm]['component_contribution'][0]['group_vector']
                except KeyError:
                    pass
                try:
                    pmap_species = self.mnxm_dG[mnxm]['component_contribution'][0]['pmap']['species']
                except KeyError:
                    raise KeyError
                #if compound_index==None and group_vector==None and pmap_species==None:
                    #exit the component_conribution condition if all are none
                #    continue
                return compound_index, group_vector, pmap_species
            #if you cannot find the component in the component_contribution then select the alberty dataset
            elif 'alberty' in self.mnxm_dG[mnxm] and self.mnxm_dG[mnxm]['alberty']:
                if len(self.mnxm_dG[mnxm]['alberty'])==1:
                    return None, None, self.mnxm_dG[mnxm]['alberty'][0]['species']
                else:
                    return False
            else:
                logging.warning('There are no valid dictionnary of precalculated dG for '+str(mnxm))
                raise KeyError
        else:
            logging.warning('There are no '+str(mnxm)+' in self.mnxm_dG')
            raise KeyError


    #######################################################
    ################ PUBLIC FUNCTIONS #####################
    #######################################################

    ################# ON THE FLY #######################

    #must make sure that there is no KEGG id and that there is a SMILES description
    #TODO export the matrices
    #TODO add the option of accepting InChI strings as well as SMILES
    def scrt_dG0(self, srct_type, srct_string, stoichio):
        """ Decompose a SMILES string
        #Warning -- the dimensions of x and g are not the same as the compound_to_matrix function
        calculate pKas of the target and intermediates using cxcalc
        """
        molecule = None
        if srct_type=='smiles':
            molecule = pybel.readstring('smiles', srct_string)
        elif srct_type=='inchi':
            molecule = pybel.readstring('inchi', srct_string)
        else:
            logging.error('Must input a valid molecular structure string')
            raise LookupError
        inchi = molecule.write('inchi').strip()
        inchi_key = molecule.write("inchikey").strip()
        if not inchi_key:
            logging.error('Molecule with no explicit structure: '+str(srct_string))
            raise LookupError
        #compute pKas
        try:
            p_kas, major_ms_smiles = component_contribution.chemaxon.get_dissociation_constants(inchi)
        except:
            logging.error('ChemAxon has encountered an error')
            raise LookupError
        p_kas = sorted([pka for pka in p_kas if self.min_pH<pka<self.max_pH], reverse=True)
        molecule = pybel.readstring('smi', major_ms_smiles)
        atom_bag, major_ms_charge = component_contribution.compound.atom_bag_and_charge(molecule)
        n_o_p = atom_bag.get('H', 0)
        n_species = len(p_kas)+1
        if not p_kas:
            major_microspecies = 0
        else:
            major_microspecies = len([1 for pka in p_kas if pka>7])
        number_of_protons = []
        charges = []
        for i in range(n_species):
            charges.append((i-major_microspecies)+major_ms_charge)
            number_of_protons.append((i-major_microspecies)+n_o_p)
        try:
            g = self.decomposer.smiles_to_groupvec(molecule.write('smiles')).as_array()
        except component_contribution.inchi2gv.GroupDecompositionError as gde:
            logging.error('Cannot decompose SMILES: '+str(srct_string))
            #return None, None, None
            raise LookupError
        ### using equilibrator training data
        X = np.array(np.zeros((self.cc_preprocess['C1'].shape[0], 1)), ndmin=2)
        G = np.array(np.zeros((self.cc_preprocess['C3'].shape[0], 1)), ndmin=2)
        ############################### calculate dG0_r
        G[:len(self.group_names), 0] = g #WARNING: inspired from component_contribution, here using equilibrator data
        dG0_cc = X.T @ self.cc_preprocess['v_r'] + \
                 G.T @ self.cc_preprocess['v_g']
        dG0_cc = dG0_cc[0]
        #### _dG0_vector
        dG0s = -np.cumsum([0] + p_kas) * self.R * self.temperature * np.log(10)
        # dG0' = dG0 + nH * (R T ln(10) pH + DH) - charge^2 * DH
        pseudoisomers = np.vstack([dG0s, np.array(number_of_protons), np.array(charges)]).T
        dG0_vector = pseudoisomers[:, 0] + \
            pseudoisomers[:, 1]*(self.R*self.temperature*np.log(10)*self.pH+self.debye_huckel) - \
            pseudoisomers[:, 2]**2 * self.debye_huckel
        #### _transform
        trans = -self.R * self.temperature * logsumexp(dG0_vector / (-self.R * self.temperature))
        #### _ddG
        ddG = None
        if 0==major_microspecies:
            ddG = 0
        elif 0<self.I:
            ddG = sum(p_kas[0:major_microspecies]) * self.R * self.temperature * np.log(10)
        else:
            ddG = -sum(p_kas[major_microspecies:0]) * self.R * self.temperature * np.log(10)
        ddG0_forward = trans+ddG
        dG0 = dG0_cc+ddG0_forward
        toRet = stoichio*dG0
        if type(toRet)==np.ndarray:
            toRet = toRet[0].astype(float)
        return toRet, X, G, {'atom_bag': atom_bag, 'p_kas': p_kas, 'major_ms': major_microspecies, 'number_of_protons': number_of_protons, 'charges': charges}


    ################### PRECALCULATED ########################


    def compound_dG0(self, species_list, compound_index, group_vector, stoichio):
        """Get a detla-deltaG estimate for this group of species.
            i.e., this is the difference between the dG0 and the dG'0, which
            only depends on the pKa of the pseudoisomers, but not on their
            formation energies.
        """
        # Compute per-species transforms, scaled down by R*T.
        dG0_vec = []
        for cmp_spe in species_list:
            try:
                #Transform this individual estimate to difference conditions.
                sqrt_I = np.sqrt(self.I)
                ddG_prime = 0
                # add the potential related to the pH
                if cmp_spe['nH']>0:
                    ddG_prime += cmp_spe['nH']*self.RTlog10*self.pH
                # add the potential related to the ionic strength
                ddG_prime -= self.debye_hueckle_a*(cmp_spe['z']**2-cmp_spe['nH'])*sqrt_I/(1.0+self.debye_hueckle_b*sqrt_I)
                # add the potential related to the Mg ions
                if cmp_spe['nMg']>0:
                    ddG_prime += cmp_spe['nMg']*(self.RTlog10*self.pMg-self.mg_formation_energy)
                dG0_vec.append(cmp_spe['dG0_f']+ddG_prime)
            except KeyError:
                #TODO: choose wether to continue the calculation although some of the species_list is missing parameters
                #I think continue is fine ---> check against equilibrator
                logging.warning('Some species_list parameters are missing')
                continue
                #raise KeyError
        dG0_vec = np.array(dG0_vec)
        # Numerical issues: taking a sum of exp(v) for |v| quite large.
        # Use the fact that we take a log later to offset all values by a
        # constant (the minimum value).
        if len(dG0_vec)>0:
            dG0 = -self.RT*np.logaddexp.reduce((-1.0/self.RT)*dG0_vec)
        else:
            logging.warning('Compounds with unspecific structure')
            raise KeyError
        #all that is required with this return is a sum
        ###### compute X and G
        x = np.array(np.zeros((self.cc_preprocess['C1'].shape[0], 1)))
        g = np.array(np.zeros((self.cc_preprocess['C3'].shape[0], 1)))
        if not compound_index==None and not group_vector==None:
            for g_ind, g_count in group_vector:
                g[g_ind, 0] += g_count
            x[compound_index, 0] = 1
        #return dG0, x, g
        return dG0*stoichio, x*stoichio, g*stoichio

    #i_to --> major_ms
    #i_from --> 0
    #TODO: make this work--- right now too many errors
    def compound_dG0_prime(self, dG0, p_kas, number_of_protons, charges, stoichio, i_to, i_from=0):
        ## _dG0_prime_vector ## 
        #Calculates the difference in kJ/mol between dG'0 and dG0 of the MS with the least hydrogens (dG0[0])
        #if not self.inchi:
        #    return 0
        if not p_kas:
            dG0s = np.zeros((1, 1))
        else:
            dG0s = -np.cumsum([0]+p_kas)*self.R*self.temperature*np.log(10)
            dG0s = dG0s
        pseudoisomers = np.vstack([dG0s, np.array(number_of_protons), np.array(charges)]).T
        dG0_prime_vector = pseudoisomers[:, 0] + \
            pseudoisomers[:, 1]*(self.R*self.temperature*np.log(10)*self.pH+self.debye_huckel) - \
            pseudoisomers[:, 2]**2*self.debye_huckel
        ## _ddG ##
        #Calculate the difference in kJ/mol between two MSs
        ddG = 0
        if not (0 <= i_from <= len(p_kas)):
            raise ValueError('MS index is out of bounds: 0 <= %d <= %d' % (
                i_from, len(p_kas)))
        if not (0 <= i_to <= len(p_kas)):
            raise ValueError('MS index is out of bounds: 0 <= %d <= %d' % (
                i_to, len(p_kas)))
        if i_from == i_to:
            ddG = 0
        elif i_from < i_to:
            ddG = sum(p_kas[i_from:i_to])*self.R*self.temperature*np.log(10)
        else:
            ddG = -sum(p_kas[i_to:i_from])*self.R*self.temperature*np.log(10)
        ##### transform
        #-self.R*self.temperature*logsumexp(dG0_prime_vector/(-self.R*self.temperature))
        #return dG0_prime_vector
        dG0_prime_vector = -self.R*self.temperature*logsumexp(dG0_prime_vector/(-self.R*self.temperature)) 
        return dG0+stoichio*(dG0_prime_vector+ddG)


    def sigma_matrix(self, X, G):
        """Calculate the uncertainty  of a precalculated compound using matrices X and G
        """
        #Calculate the uncertainty of a precalculated compound or group of compounds 
        U = X.T @ self.cc_preprocess['C1'] @ X + \
            X.T @ self.cc_preprocess['C2'] @ G + \
            G.T @ self.cc_preprocess['C2'].T @ X + \
            G.T @ self.cc_preprocess['C3'] @ G
        #Below is terrible score, if we want to calculate the dG0 using the matrix
        #dG0_cc = X.T @ self.cc_preprocess['v_r'] + G.T @ self.cc_preprocess['v_g']
        return np.sqrt(U[0, 0])#, dG0_cc


    #### select either precalculated or structure dG for the calculation of dG for a pathway ####

    #TODO: change to update the JSON object instead
    def rp_paths_dG0(self, rp_paths, rp_smiles):
        """Return the paths dG from precalculated and on-the-fly
        """
        #TODO: save the caluclated MNXM and save them to the cache to speed up future calculations
        for path_id in rp_paths:
            path_dG = 0.0
            path_count = 0
            path_num_reactions = sum([len(rp_paths[path_id]['path'][step_id]['steps']) for step_id in rp_paths[path_id]['path']])
            X_path = np.array(np.zeros((self.cc_preprocess['C1'].shape[0], path_num_reactions)), ndmin=2)
            G_path = np.array(np.zeros((self.cc_preprocess['C3'].shape[0], path_num_reactions)), ndmin=2)
            for step_id in rp_paths[path_id]['path']:
                reaction_dG = 0.0
                react_count = 0
                X_reaction = np.array(np.zeros((self.cc_preprocess['C1'].shape[0], len(rp_paths[path_id]['path'][step_id]['steps']))), ndmin=2)
                G_reaction = np.array(np.zeros((self.cc_preprocess['C3'].shape[0], len(rp_paths[path_id]['path'][step_id]['steps']))), ndmin=2)
                for cmp_id in rp_paths[path_id]['path'][step_id]['steps']:
                    X = None
                    G = None
                    compound_dG = None
                    #check that the deprecated mnxm is not used
                    try:
                        mnxm = self.deprecatedMNXM_mnxm[cmp_id]
                    except KeyError:
                        mnxm = cmp_id
                    try:
                        logging.info('######### '+str(mnxm)+' ##############')
                        compound_index, group_vector, species_list = self._select_mnxm_dG(mnxm)
                        compound_dG, X, G = self.compound_dG0(species_list,
                                                              compound_index,
                                                              group_vector,
                                                              rp_paths[path_id]['path'][step_id]['steps'][cmp_id]['stoichiometry'])
                    except KeyError:
                        logging.warning(str(cmp_id)+' cannot be found in the precalculated mnxm_dG ('+str(path_id)+', '+str(step_id)+', '+str(cmp_id)+')')
                        try:
                            if cmp_id in self.calculated_dG:
                                compound_dG = self.calculated_dG[cmp_id]['compound_dG']
                                X = self.calculated_dG[cmp_id]['X']
                                G = self.calculated_dG[cmp_id]['G']
                                cmp_info = self.calculated_dG[cmp_id]['cmp_info']
                            else:
                                compound_dG, X, G, cmp_info = self.scrt_dG0('smiles', rp_smiles[cmp_id], rp_paths[path_id]['path'][step_id]['steps'][cmp_id]['stoichiometry'])
                                #save it for the next iteractions to avoid duplicates
                                self.calculated_dG[cmp_id] = {'compound_dG': compound_dG, 'X': X, 'G': G, 'cmp_info': cmp_info}
                        except KeyError:
                            logging.warning(str(cmp_id)+' cannot be found in the RP2paths compounds.txt ('+str(path_id)+', '+str(step_id)+', '+str(cmp_id)+')')
                        except LookupError:
                            logging.warning(str(cmp_id)+' cannot calculate dG using component contribution ('+str(path_id)+', '+str(step_id)+', '+str(cmp_id)+')')
                    if type(X)==np.ndarray and type(G)==np.ndarray:
                        X_reaction[:, react_count:react_count+1] = X
                        G_reaction[:, react_count:react_count+1] = G
                        X_path[:, path_count:path_count+1] = X
                        G_path[:, path_count:path_count+1] = G
                        rp_paths[path_id]['path'][step_id]['steps'][cmp_id]['dG_uncertainty'] = self.sigma_matrix(X, G)
                    if not compound_dG==None:
                        rp_paths[path_id]['path'][step_id]['steps'][cmp_id]['dG'] = compound_dG
                        reaction_dG += compound_dG
                        path_dG += compound_dG
                    react_count += 1
                    path_count += 1
                    logging.info('############################')
                rp_paths[path_id]['path'][step_id]['dG'] = reaction_dG
                rp_paths[path_id]['path'][step_id]['dG_uncertainty'] = self.sigma_matrix(X_reaction, G_reaction)
            rp_paths[path_id]['dG'] = path_dG
            rp_paths[path_id]['dG_uncertainty'] = self.sigma_matrix(X_path, G_path)



    #TODO: save the calculated InChI dG and uncertainty to be used once again 
    def rp_paths_dG0_dGPrime(self, rp_paths, rp_smiles):
        """Return the paths dG from precalculated and on-the-fly
        """
        #TODO: save the caluclated MNXM and save them to the cache to speed up future calculations
        for path_id in rp_paths:
            path_dG = 0.0
            path_dG_prime = 0.0
            path_count = 0
            path_num_reactions = sum([len(rp_paths[path_id]['path'][step_id]['steps']) for step_id in rp_paths[path_id]['path']])
            X_path = np.array(np.zeros((self.cc_preprocess['C1'].shape[0], path_num_reactions)), ndmin=2)
            G_path = np.array(np.zeros((self.cc_preprocess['C3'].shape[0], path_num_reactions)), ndmin=2)
            for step_id in rp_paths[path_id]['path']:
                reaction_dG = 0.0
                reaction_dG_prime = 0.0
                react_count = 0
                X_reaction = np.array(np.zeros((self.cc_preprocess['C1'].shape[0], 
                    len(rp_paths[path_id]['path'][step_id]['steps']))), ndmin=2)
                G_reaction = np.array(np.zeros((self.cc_preprocess['C3'].shape[0], 
                    len(rp_paths[path_id]['path'][step_id]['steps']))), ndmin=2)
                for cmp_id in rp_paths[path_id]['path'][step_id]['steps']:
                    X = None
                    G = None
                    cmp_dG = None
                    cmp_dG_prime = None
                    #check that the deprecated mnxm is not used
                    try:
                        mnxm = self.deprecatedMNXM_mnxm[cmp_id]
                    except KeyError:
                        mnxm = cmp_id
                    try:
                        logging.info('######### '+str(mnxm)+' ##############')
                        compound_index, group_vector, species_list = self._select_mnxm_dG(mnxm)
                        logging.info('compound_index: '+str(compound_index))
                        logging.info('group_vector: '+str(group_vector))
                        logging.info('species_list: '+str(species_list))
                        cmp_dG, X, G = self.compound_dG0(species_list,
                                        compound_index,
                                        group_vector,
                                        rp_paths[path_id]['path'][step_id]['steps'][cmp_id]['stoichiometry'])
                        try:
                            cmp_dG_prime = self.compound_dG0_prime(cmp_dG, 
                                            rp_paths[path_id]['path'][step_id]['steps'][cmp_id]['p_kas'], 
                                            rp_paths[path_id]['path'][step_id]['steps'][cmp_id]['number_of_protons'],
                                            rp_paths[path_id]['path'][step_id]['steps'][cmp_id]['charges'],
                                            rp_paths[path_id]['path'][step_id]['steps'][cmp_id]['stoichiometry'],
                                            rp_paths[path_id]['path'][step_id]['steps'][cmp_id]['major_ms'])
                        except KeyError:
                            pass
                    except KeyError:
                        logging.warning(str(cmp_id)+' cannot be found in the precalculated mnxm_dG ('+str(path_id)+', '+str(step_id)+', '+str(cmp_id)+')')
                        try:
                            cmp_dG, X, G, cmp_info = self.scrt_dG0('smiles', 
                                                rp_smiles[cmp_id], 
                                                rp_paths[path_id]['path'][step_id]['steps'][cmp_id]['stoichiometry'])
                            cmp_dG_prime = self.compound_dG0_prime(cmp_dG, 
                                                cmp_info['p_kas'], 
                                                cmp_info['number_of_protons'],
                                                cmp_info['charges'],
                                                rp_paths[path_id]['path'][step_id]['steps'][cmp_id]['stoichiometry'],
                                                cmp_info['major_ms'])
                        except KeyError:
                            logging.warning(str(cmp_id)+' cannot be found in the RP2paths compounds.txt ('+str(path_id)+', '+str(step_id)+', '+str(cmp_id)+')')
                        except LookupError:
                            logging.warning(str(cmp_id)+' cannot calculate dG using component contribution ('+str(path_id)+', '+str(step_id)+', '+str(cmp_id)+')')
                    if type(X)==np.ndarray and type(G)==np.ndarray:
                        X_reaction[:, react_count:react_count+1] = X
                        G_reaction[:, react_count:react_count+1] = G
                        X_path[:, path_count:path_count+1] = X
                        G_path[:, path_count:path_count+1] = G
                        rp_paths[path_id]['path'][step_id]['steps'][cmp_id]['dG_uncertainty'] = self.sigma_matrix(X, G)
                    if not cmp_dG==None and not cmp_dG_prime==None:
                        rp_paths[path_id]['path'][step_id]['steps'][cmp_id]['dG'] = cmp_dG
                        reaction_dG += cmp_dG
                        reaction_dG_prime += cmp_dG_prime
                        path_dG += cmp_dG
                        path_dG_prime += cmp_dG_prime
                    react_count += 1
                    path_count += 1
                    logging.info('############################')
                rp_paths[path_id]['path'][step_id]['dG'] = reaction_dG
                rp_paths[path_id]['path'][step_id]['dG_prime'] = reaction_dG_prime
                rp_paths[path_id]['path'][step_id]['dG_uncertainty'] = self.sigma_matrix(X_reaction, G_reaction)
            rp_paths[path_id]['dG'] = path_dG
            rp_paths[path_id]['dG_prime'] = path_dG_prime
            rp_paths[path_id]['dG_uncertainty'] = self.sigma_matrix(X_path, G_path)


    def isBalanced():
        """Function borrowed from the component contribution that checks is the per-atom
        difference in a reaction is balanced
        """
        #TODO: check the difference between equilibrator and component_contribution
        return False

    def reversibilityIndex():
        """Quantitative measure for the reversibility of a reaction by taking into consideration the concentration of the substrate and products
        """
        return False

    """
        IDEA: run a small scale kinetic model of the heterologous pathway by estimating from the FBA analysis the concentration of the sink molecules. We can estimate the steady state concentration and then calculate other features such as the reversibility index of the reactions within the pathway.
        Would this not be the same as deriving the concentration from the flux
    """

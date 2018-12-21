import numpy as np
import logging
import pybel
import json
from scipy.special import logsumexp
import copy

#local package
import component_contribution


class Thermodynamics:
    """Combination of equilibrator and group_contribution analysis to caclualte the thermodymaics of the individual 
    metabolic pathways output from RP2paths and thus RetroPath2.0
    """
    def __init__(self, mnxm_dG, cc_preprocess, pH=7.0, pMg=14.0, I=0.1, temperature=298.15):
        #self.compound_dict = compound_dict
        self.mnxm_dG = mnxm_dG
        self.cc_preprocess = cc_preprocess
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


    ################ PRIVATE FUNCTIONS ###################


    def _select_mnxm_dG(mnxm):
        """Given that there can be multiple precalculated dG for a given molecule (MNXM)
        this function gives the priority to a particular order for a given MNXM ID
        alberty>smallest(KEGG_ID)>other(KEGG_ID)
        """
        if mnxm in self.mnxm_dG:
            if 'alberty' in self.mnxm_dG[mnxm] and self.mnxm_dG[mnxm]['alberty']:
                if len(self.mnxm_dG[mnxm]['alberty'])==1:
                    return None, None, self.mnxm_dG[mnxm]['alberty'][0]['species']
                else:
                    return False
            elif 'component_contribution' in self.mnxm_dG[mnxm] and self.mnxm_dG[mnxm]['component_contribution']:
                if len(self.mnxm_dG[mnxm]['component_contribution'])==1:
                    return self.mnxm_dG[mnxm]['component_contribution'][0]['compound_index'], self.mnxm_dG[mnxm]['component_contribution'][0]['group_vector'], self.mnxm_dG[mnxm]['component_contribution'][0]['pmap']['species']
                elif:
                    #select the lowest CID and make sure that there
                    toRet = {'CID': 'C99999'}
                    for cmp_dict in self.mnxm_dG[mnxm]['component_contribution']:
                        if int(cmp_dict['CID'][1:])<int(toRet['CID'][1:]):
                            toRet = cmp_dict
                    return toRet['compound_index'], toRet['group_vector'], toRet['pmap']['species']
            else:
                logging.warning('There are no valid dictionnary of precalculated dG for '+str(mnxm))
                raise KeyError
        else:
            logging.warning('There are no '+str(mnxm)+' in self.mnxm_dG')
            raise KeyError


    ################ PUBLIC FUNCTIONS #####################


    ################# MOLECULAR STRUCTURE #######################


    #must make sure that there is no KEGG id and that there is a SMILES description
    #TODO export the matrices
    #TODO generate the uncertainty here
    #TODO add the option of accepting InChI strings as well as SMILES
    def scrt_dG0(self, srct_type, srct_string, stochio):
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
        #### _dG0_prime_vector
        dG0s = -np.cumsum([0] + p_kas) * self.R * self.temperature * np.log(10)
        # dG0' = dG0 + nH * (R T ln(10) pH + DH) - charge^2 * DH
        pseudoisomers = np.vstack([dG0s, np.array(number_of_protons), np.array(charges)]).T
        dG0_prime_vector = pseudoisomers[:, 0] + \
            pseudoisomers[:, 1]*(self.R*self.temperature*np.log(10)*self.pH+self.debye_huckel) - \
            pseudoisomers[:, 2]**2 * self.debye_huckel
        #### _transform
        trans = -self.R * self.temperature * logsumexp(dG0_prime_vector / (-self.R * self.temperature))
        #### _ddG
        ddG = None
        if 0 == major_microspecies:
            ddG = 0
        elif 0 < self.I:
            ddG = sum(p_kas[0:major_microspecies]) * self.R * self.temperature * np.log(10)
        else:
            ddG = -sum(p_kas[major_microspecies:0]) * self.R * self.temperature * np.log(10)
        ddG0_forward = trans+ddG
        dG0_prime = dG0_cc+ddG0_forward
        toRet = stochio*dG0_prime
        if type(toRet)==np.ndarray:
            toRet = toRet[0].astype(float)
        return toRet, X, G


    ################### PRECALCULATED ########################


    def compound_dG0(self, species_list, compound_index, group_vector, stochio):
        """Get a detla-deltaG estimate for this group of species.
            i.e., this is the difference between the dG0 and the dG'0, which
            only depends on the pKa of the pseudoisomers, but not on their
            formation energies.          
        """
        # Compute per-species transforms, scaled down by R*T.
        dG0_prime_vec = []
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
                dG0_prime_vec.append(cmp_spe['dG0_f']+ddG_prime)
            except KeyError:
                #TODO: pass the mnxm to print where is the error
                logging.warning('Some species_list parameters are missing')
                continue
        dG0_prime_vec = np.array(dG0_prime_vec)
        # Numerical issues: taking a sum of exp(v) for |v| quite large.
        # Use the fact that we take a log later to offset all values by a
        # constant (the minimum value).
        if len(dG0_prime_vec)>0:
            dG0_f_prime = -self.RT*np.logaddexp.reduce((-1.0/self.RT)*dG0_prime_vec)
        else:
            logging.warning('Compounds with unspecific structure')
            raise KeyError
        #all that is required with this return is a sum
        ###### compute X and G
        x = np.array(np.zeros((self.cc_preprocess['C1'].shape[0], 1)))
        g = np.array(np.zeros((self.cc_preprocess['C3'].shape[0], 1)))
        for g_ind, g_count in group_vector:
            g[g_ind, 0] += g_count
        x[compound_index, 0] = 1
        return dG0_f_prime*stochio, stochio*x, stochio*g

    
    def sigma_matrix(self, X, G)
        #Calculate the uncertainty of a precalculated compound or group of compounds 
        U = X.T @ self.cc_preprocess['C1'] @ X + \
            X.T @ self.cc_preprocess['C2'] @ G + \
            G.T @ self.cc_preprocess['C2'].T @ X + \
            G.T @ self.cc_preprocess['C3'] @ G 
        #Below is terrible score, if we want to calculate the dG0 using the matrix
        #dG0_cc = X.T @ self.cc_preprocess['v_r'] + G.T @ self.cc_preprocess['v_g']
        return np.sqrt(U[0, 0])#, dG0_cc
    

    #### select either precalculated or structure dG for the calculation of dG for a pathway ####

    
    def rp_paths_dG0(self, input_rp_paths, rp_smiles):
        """Return the paths dG from precalculated and on-the-fly 
        """
        rp_paths = copy.deepcopy(input_rp_paths)
        for path in rp_paths:
            path_dG = 0.0
            path_count = 0
            path_num_reactions = sum([len(react['right'])+len(react['left']) for react in path])
            X_path = np.array(np.zeros((self.cc_preprocess['C1'].shape[0], path_num_reactions)), ndmin=2)
            G_path = np.array(np.zeros((self.cc_preprocess['C3'].shape[0], path_num_reactions)), ndmin=2)
            for react in path:
                reaction_dG = 0.0
                react_count = 0
                X_reaction = np.array(np.zeros((self.cc_preprocess['C1'].shape[0], 
                                                len(react['right'])+len(react['left']))), ndmin=2) 
                G_reaction = np.array(np.zeros((self.cc_preprocess['C3'].shape[0], 
                                                len(react['right'])+len(react['left']))), ndmin=2)
                for lr, stochio_mult in zip(['left', 'right'], [-1, 1]):
                    for mnxm in react[lr]:
                        try:
                            compound_index, group_vector, species_list = self._select_mnxm_dG(mnxm)
                            compound_dG, X, G = self.compound_dG0(species_list, 
                                                                  compound_index, 
                                                                  group_vector, 
                                                                  react[lr][mnxm]*stochio_mult)
                        except KeyError:
                            try:
                                compound_dG, X, G = scrt_dG0('smiles', rp_smiles[mnxm], -react['left'][mnxm])
                            except KeyError:
                                logging.warning('This ID ('+str(l)+') cannot be found in the RP2paths compounds.txt')
                            except LookupError:
                                logging.warning('Cannot calculate component contribution for ID ('+str(l)+')')
                        if not X==None and not G==None:
                            X_reaction[react_count] = X
                            G_reaction[react_count] = G
                            X_path[path_count] = X
                            G_path[path_count] = G
                        reaction_dG += compound_dG
                        react_count += 1
                        path_count += 1
                react['dG'] = compound_dG
                react['dG_uncertainty'] = sigma_matrix(X, G)
                path_dG += react_dG            
            path.append({'dG': path_dG, 'dG_uncertainty': sigma_matrix(X_path, G_path)})
        return rp_paths






























    def rp_paths_dG0(self, input_rp_paths):
        """return the paths dG from precalculated and 
        """
        rp_paths = copy.deepcopy(input_rp_paths)
        smiles_calculated = {}
        inchi_calculated = {}
        path_count = 0
        for path_id in rp_paths:
            path_dG = 0.0
            all_step_path = sum([len(rp_paths[path_id]['path'][i]['reaction']) for i in rp_paths[path_id]['path']])
            X_path = np.array(np.zeros((self.cc_preprocess['C1'].shape[0], all_step_path)), ndmin=2)
            G_path = np.array(np.zeros((self.cc_preprocess['C3'].shape[0], all_step_path)), ndmin=2)
            for react_count, rule_id in enumerate(rp_paths[path_id]['path']):
                reaction_dG = 0.0
                X_reaction = np.array(np.zeros((self.cc_preprocess['C1'].shape[0], 
                                                len(rp_paths[path_id]['path'][rule_id]))), ndmin=2)
                G_reaction = np.array(np.zeros((self.cc_preprocess['C3'].shape[0], 
                                                len(rp_paths[path_id]['path'][rule_id]))), ndmin=2)
                for compound in rp_paths[path_id]['path'][rule_id]['reaction']:
                    try:
                        cmp_dG = self._select_mnxm_dG()
                    except KeyError:
                        



                    ####### USE THE KEGG ID ######
                    try:
                        cmp_dG = self.compound_dG0(compound['kegg'], compound['stochio'])
                        X, G = self.compound_vectors(compound['kegg'], compound['stochio'])
                        cmp_std = self.matrix_uncertainty(X, G)
                    except KeyError:
                        try:
                            ###### USE SMILES ###########
                            if compound['smiles'] in smiles_calculated:
                                cmp_dG = smiles_calculated[compound['smiles']]['cmp_dG']
                                X = smiles_calculated[compound['smiles']]['X']
                                G = smiles_calculated[compound['smiles']]['G']
                                cmp_std = self.matrix_uncertainty(X, G) 
                            else:
                                cmp_dG, X, G = self.dG0_prime('smiles', compound['smiles'], compound['stochio'])
                                cmp_std = self.matrix_uncertainty(X, G) 
                                smiles_calculated[compound['smiles']] = {'cmp_dG': cmp_dG, 'X': X, 'G': G}
                        except (KeyError, LookupError) as e:
                            try:
                                ##### USE INCHI ###########
                                if compound['inchi'] in inchi_calculated:
                                    cmp_dG = inchi_calculated[compound['inchi']]['cmp_dG']
                                    X = inchi_calculated[compound['inchi']]['X']
                                    G = inchi_calculated[compound['inchi']]['G']
                                    cmp_std = self.matrix_uncertainty(X, G) 
                                else:
                                    cmp_dG, X, G = self.dG0_prime('inchi',compound['inchi'],compound['stochio']) 
                                    cmp_std = self.matrix_uncertainty(X, G) 
                                    inchi_calculated[compound['inchi']] = {'cmp_dG': cmp_dG, 'X': X, 'G': G}
                            except (KeyError, LookupError) as e:
                                # if all fails then loop through the next compound adding nothing
                                path_count += 1
                                continue
                    #for the compound
                    compound['dG'] = cmp_dG
                    compound['dG_std'] = cmp_std
                    #for the reaction
                    reaction_dG += cmp_dG
                    X_reaction[:, react_count:react_count+1] = X
                    G_reaction[:, react_count:react_count+1] = G
                    #for the path
                    path_dG += cmp_dG
                    X_path[:, path_count:path_count+1] = X
                    G_path[:, path_count:path_count+1] = G
                    path_count += 1
                rp_paths[path_id]['path'][rule_id]['dG'] = reaction_dG
                rp_paths[path_id]['path'][rule_id]['dG_std'] = self.matrix_uncertainty(X_reaction, G_reaction)
            rp_paths[path_id]['dG'] = path_dG
            rp_paths[path_id]['dG_std'] = self.matrix_uncertainty(X_path, G_path)
        return rp_paths 






















    ########################### WORK IN PROGRESS #####################

    ################## Compound to matrix and calculate uncertainty ################
    """
    #Calulate the matrices to calculate the uncertainty parameter
    def compound_vectors(self, kegg_id, stochio):
        ###get_stoich_vector
        x = np.array(np.zeros((self.cc_preprocess['C1'].shape[0], 1)))
        try:
            compound_index = self.compound_dict[kegg_id]['compound_index']
        except KeyError:
            logging.error('This kegg_id ('+str(kegg_id)+') does not seem to exist in self.compound_dict or does not contain a compound index')
            raise KeyError
        x[compound_index, 0] = 1
        return stochio*x

    def group_incidence_vector(self, kegg_id, stochio):
        g = np.array(np.zeros((self.cc_preprocess['C3'].shape[0], 1)))
        try:
            group_vector = self.compound_dict[kegg_id]['group_vector']
        except KeyError:
            logging.error('This kegg_id ('+str(kegg_id)+') does not seem to exist or does not contain group_vector')
            raise KeyError
        for g_ind, g_count in group_vector:
            g[g_ind, 0] += g_count
        return stochio*g 

    def 

    def reaction_to_matrix(self, reaction):
        x = np.array(np.zeros((self.cc_preprocess['C1'].shape[0], 1)))
        g = np.array(np.zeros((self.cc_preprocess['C3'].shape[0], 1)))
        for compound in reaction:
            if compound['kegg']:
                x += self.compound_vectors(compound['kegg'], compound['stochio'])
                g += self.group_incidence_vector(compound['kegg'], compound['stochio'])
            else:
                logging.warning('Passing an empty kegg to the reaction_to_matrix')
        ## not sure about this transformation
        X = np.array(np.zeros((self.cc_preprocess['C1'].shape[0], 1)), ndmin=2)
        G = np.array(np.zeros((self.cc_preprocess['C3'].shape[0], 1)), ndmin=2)
        X[:, 0:1] = x
        G[:, 0:1] = g
        return X, G
    """


    """
    def reaction_uncertainty(self, reaction):
        X,G = self.reaction_to_matrix(reaction)
        return self.matrix_uncertainty(X, G)


    def path_uncertainty(self, path):
        for reaction in path:
            return False
    """

    ########################### REACTIONS #####################

    ###uselesss since we are passing

    #for standard 1M concentrations
    #TODO return the uncertainty and error margin
    #TODO: see if you can generate the same results but using the matrices instead of the stored values
    #   -> the idea is to cobine the matrices for a reaction or set of reactions
    def reaction_dG0_prime(self, reaction):
        #1)test the validity of the input
        dG0_prime = 0.0
        for compound in reaction:
            if compound['kegg'] and not compound['smiles']:
                if compound['kegg'] in self.cc_compounds:
                    #TODO: change this to compound or change the below to smiles
                    dG0_prime += self.compound_dG0(compound['kegg'], compound['stochio'])
                else:
                    logger.error('Compound '+str(compound['kegg']+' is not contained in the equilibrator list'))
            if not compound['kegg'] and compound['smiles']:
                dG0_prime += self.dG0_prime(compound['smiles'], compound['stochio'])[0]
            else:
                logging.error('The compound contains no KEGG and no SMILES info')
                return -1.0
        return dG0_prime, reaction_uncertainty(reaction)
   

    #Do we really need this?
    #TODO Change to accomodate for the different input
    def reaction_dG_prime_conc(self, reaction):
        """Based on the known concentrations of the compounds in a reaction calculate the dG
        """
        #calculate the reaction dG0_prime
        dG0_prime = self.dG0_prime(reaction)
        dG_correction = 0.0
        for compound in reaction:
            dG_correction += compound['stochio']*np.log(compound['conc'])
        return dG0_prime+(self.RT*dG_correction), reaction_uncertainty(reaction)


    #TODO add for the condition of 1mM
    def reaction_dGm_prime(self, reaction):
        """Calculate the free energy of a reaction for a 1mM and not 1M condition
        """
        dG0_prime = self.dG0_prime(reaction)
        dGm_correction = 0.0
        for compound in reaction:
            dGm_correction += self.RT*compound['stochio']*np.log(1e-3)
        return dG0_prime+dGm_correction, reaction_uncertainty(reaction)





    #TODO
    #check that there is the minimal information required for the reation to run smoothly
    #check that the KEGG id exists in the preprocessed data from equilibrium 
    def check_inputReaction_integrity(self, inputReaction):
        return False
 

    #TODO
    def check_reaction_is_balanced(self, inputReaction):
        return False


    #TODO
    def atom_bag(self, compound):
        return False


    def reaction_atom_balance(self, reaction):
        return False

    def calculate_reversibility_index_from_dG0_prime(self):
        return False

    def E0_prime(self):
        return False

    def check_reaction_balance(self):
        #get_atom_bag
        atom_bag = {}

        #elements, Ematrix
        #get_element_matrix

    


    ############################## Matrices ###################################

    def scrt_vector(self, srct_type, srct_string, stochio):
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
            return None, None
        inchi = molecule.write('inchi').strip()
        inchi_key = molecule.write("inchikey").strip()
        if not inchi_key:
            logging.error('Molecule with no explicit structure: '+str(srct_string))
            return None, None
        #compute pKas
        try:
            p_kas, major_ms_smiles = component_contribution.chemaxon.get_dissociation_constants(inchi)
        except:
            logging.error('ChemAxon has encountered an error')
            #logging.error(e)
            return None, None
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
            return None, None
        ### using equilibrator training data
        x = np.array(np.zeros((self.cc_preprocess['C1'].shape[0], 1)), ndmin=2)
        ###### calculate dG0_r
        return stochio*x, stochio*g

    """
    #Calulate the matrices to calculate the uncertainty parameter
    def compound_vectors(self, kegg_id, stochio):
        ###get_stoich_vector
        x = np.array(np.zeros((self.cc_preprocess['C1'].shape[0], 1)))
        g = np.array(np.zeros((self.cc_preprocess['C3'].shape[0], 1)))
        try:
            compound_index = self.compound_dict[kegg_id]['compound_index']
        except KeyError:
            logging.warning('This kegg_id ('+str(kegg_id)+') does not seem to exist in self.compound_dict or does not contain a compound index')
            return None, None
        try:
            for g_ind, g_count in self.compound_dict[kegg_id]['group_vector']:
                g[g_ind, 0] += g_count
        except KeyError:
            logging.warning('The self.compound_dict does not contain group_vector for '+str(kegg_id)+' or does not exist')
            return None, None
        x[compound_index, 0] = 1
        return stochio*x, stochio*g
    """

    def reaction_vector(self, reaction):
        x_r = np.array(np.zeros((self.cc_preprocess['C1'].shape[0], 1)), ndmin=2)
        g_r = np.array(np.zeros((self.cc_preprocess['C3'].shape[0], 1)), ndmin=2)
        #x_r = np.array(np.zeros((self.cc_preprocess['C1'].shape[0], 1)))
        #g_r = np.array(np.zeros((self.cc_preprocess['C3'].shape[0], 1)))
        for compound in reaction:
            print(compound)
            x = None
            g = None
            used = 0
            if compound['kegg']:
                x, g = self.compound_vectors(compound['kegg'], compound['stochio'])
                used = 1
            if not compound['smiles']==None and not type(x)==np.ndarray and not type(g)==np.ndarray:
                x, g = self.scrt_vector('smiles', compound['smiles'], compound['stochio'])
                used = 2
            if not compound['inchi']==None and not type(x)==np.ndarray and not type(g)==np.ndarray:
                x, g = self.scrt_vector('inchi', compound['inchi'], compound['stochio'])
                used = 3
            if not type(x)==np.ndarray and not type(g)==np.ndarray:
                logging.warning('Cannot calculate vector for '+str(compound['mnx']))
                return None, None
            #print(list(g_r))
            #print(list(g))
            print(g_r.shape)
            print(g.shape)
            print(used)
            x_r += x
            if g.shape[0]==163:
                g_r[:len(self.group_names), 0] += g
            else:
                g_r += g
        return x_r, g_r

    def path_vector(self, path):
        X = np.array(np.zeros((self.cc_preprocess['C1'].shape[0], len(path))), ndmin=2)
        G = np.array(np.zeros((self.cc_preprocess['C3'].shape[0], len(path))), ndmin=2)
        for i, path_id in enumerate(path):
            x, g = self.reaction_vector(path[path_id]['reaction'])
            if type(x)==np.ndarray and type(g)==np.ndarray:
                X[:, i:i+1] = x
                G[:, i:i+1] = g
            else:
                return None, None
        #return X, G
        dG0_cc = X.T @ self.cc_preprocess['v_r'] + G.T @ self.cc_preprocess['v_g']
        U = X.T @ self.cc_preprocess['C1'] @ X + \
                X.T @ self.cc_preprocess['C2'] @ G + \
                G.T @ self.cc_preprocess['C2'].T @ X + \
                G.T @ self.cc_preprocess['C3'] @ G
        #TO TEST: there is a strong chance that this does not work at all... 
        #compare with the original equilibrator_api
        return dG0_cc, np.sqrt(U[0, 0])
    #def dG0_vector(self, rp_paths):

    def matrix_dG(self, X, G):
        try:
            p_kas, major_ms_smiles = component_contribution.chemaxon.get_dissociation_constants(inchi)
        except:
            logging.error('ChemAxon has encountered an error')
            #logging.error(e)
            return None, None
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
            return None, None
        dG0_cc = X.T @ self.cc_preprocess['v_r'] + \
                 G.T @ self.cc_preprocess['v_g']
        U = X.T @ self.cc_preprocess['C1'] @ X + \
            X.T @ self.cc_preprocess['C2'] @ G + \
            G.T @ self.cc_preprocess['C2'].T @ X + \
            G.T @ self.cc_preprocess['C3'] @ G
        dG0_cc = dG0_cc[0]
        sigma_cc = np.sqrt(U[0, 0])
        #### _dG0_prime_vector
        dG0s = -np.cumsum([0] + p_kas) * self.R * self.temperature * np.log(10)
        # dG0' = dG0 + nH * (R T ln(10) pH + DH) - charge^2 * DH
        pseudoisomers = np.vstack([dG0s, np.array(number_of_protons), np.array(charges)]).T
        dG0_prime_vector = pseudoisomers[:, 0] + \
            pseudoisomers[:, 1]*(self.R*self.temperature*np.log(10)*self.pH+self.debye_huckel) - \
            pseudoisomers[:, 2]**2 * self.debye_huckel
        #### _transform
        trans = -self.R * self.temperature * logsumexp(dG0_prime_vector / (-self.R * self.temperature))
        #### _ddG
        ddG = None
        if 0 == major_microspecies:
            ddG = 0
        elif 0 < self.I:
            ddG = sum(p_kas[0:major_microspecies]) * self.R * self.temperature * np.log(10)
        else:
            ddG = -sum(p_kas[major_microspecies:0]) * self.R * self.temperature * np.log(10)
        #ddG0_forward = stochio*(trans+ddG)
        #dG0_prime = dG0_cc+ddG0_forward
        #return dG0_prime, sigma_cc
        ddG0_forward = trans+ddG
        dG0_prime = dG0_cc+ddG0_forward
        return stochio*dG0_prime, sigma_cc

    '''
    ##################### THIS IS TESTING SUITE ##################    
    def compare_scores(self):
        import equilibrator_api
        """
            Function that takes for input a list of KEGG ID's and generates the dG from the equilibrator
            data and generates it from group-contribution. 
        """
        results = {}
        not_calculated = []
        smiles_kegg_errors = []
        for kegg in kegg_mnx:
            try:
                smiles = mnx_smiles[kegg_mnx[kegg]]
                if smiles and not smiles=='NA' and kegg in self.compound_dict and not kegg in kegg_ignore:
                    try:
                        print('KEGG: '+str(kegg)+' - SMILES: '+str(smiles))
                        results[kegg] = {}
                        results[kegg]['compound'] = self.compound_dG0(kegg, 1)
                        results[kegg]['smiles'] = float(self.dG0_prime(smiles, 1)[0])
                    except GroupDecompositionError as gde:
                        print(gde)
                        del results[kegg]
                        not_calculated.append(kegg)
                        smiles_kegg_errors.append(kegg)
                else:
                    not_calculated.append(kegg)
            except KeyError:
                not_calculated.append(kegg)
                print('KEGG ID in compound_contribution is not in MNX: '+str(kegg))
        return results, not_calculated, smiles_kegg_errors
    '''

    '''
    def dG0_prime(self, inputReaction):
        dG0_r_prime = reaction_dG0_prime(inputReaction, compound_dict, pH, pMG, I)
        x, g = reaction_to_matrix(inputReaction, compound_dict)
        ##### this step seems to be useless since it is designed
        #for an array of reactions and not a single
        Nc = np.array(cc_preprocess['C1']).shape[0]
        Ng = np.array(cc_preprocess['C3']).shape[0]
        X = np.array(np.zeros((Nc, 1)), ndmin=2)
        G = np.array(np.zeros((Ng, 1)), ndmin=2)
        X[:, 0:1] = x
        G[:, 0:1] = g
        U = X.T @ cc_preprocess['C1'] @ X + \
            X.T @ cc_preprocess['C2'] @ G + \
            G.T @ cc_preprocess['C2'].T @ X + \
            G.T @ cc_preprocess['C3'] @ G
        dG0_uncertainty = np.sqrt(U[0, 0])
        return dG0_r_prime, dG0_uncertainty


    def dG_prime(inputReaction):
        dG0_r_prime, dG0_uncertainty = dG0_prime(inputReaction)
        #calculate the correction
        dG_correction = 0.0
        for i in inputReaction:
            if not inputReaction[i]['kegg']=='C00001':
                dG_correction += inputReaction[i]['stochio']*np.log(inputReaction[i]['conc'])
        return dG0_r_prime+(dG_correction*self.RT), dG0_uncertainty
    '''

























    '''
    def dG0_prime(self, inputReaction):
        x = np.array(np.zeros((np.array(self.cc_preprocess['C1']).shape[0], 1)), ndmin=2)
        g = np.array(np.zeros((np.array(self.cc_preprocess['C3']).shape[0], 1)), ndmin=2)
        for cmp_name  in inputReaction:
            if inputReaction[cmp_name]['kegg']:
                x_t, g_t = self.compound_to_matrix(inputReaction[cmp_name])
            elif not inputReaction[cmp_name]['kegg'] and inputReaction[cmp_name]['smiles']:
                x_t, g_t = self.smiles_to_matrix(inputReaction[cmp_name])
            else:
                logger.error('The input compound does not have a KEGG id nor a smiles')
        #This step seems useless, however the below calculation seems to require multi-dimensional matrices
        X = np.array(np.zeros((np.array(self.cc_preprocess['C1']).shape[0], 1)), ndmin=2)
        G = np.array(np.zeros((np.array(self.cc_preprocess['C3']).shape[0], 1)), ndmin=2)
        X[:, 0:1] = x
        G[:, 0:1] = g
        
        dG0_r_prime = X.T @ self.params['preprocess_v_r'] + G.T @ self.params['preprocess_v_g']
        U = X.T @ self.cc_preprocess['C1'] @ X + \
            X.T @ self.cc_preprocess['C2'] @ G + \
            G.T @ self.cc_preprocess['C2'].T @ X + \
            G.T @ self.cc_preprocess['C3'] @ G 
        dG0_uncertainty = np.sqrt(U[0, 0])
        return dG0_r_prime, dG0_uncertainty


    def dG_prime(self, inputReaction):
        dG0_r_prime, dG0_uncertainty = self.dG0_prime(inputReaction)
        #calculate the correction
        dG_correction = 0.0
        for i in inputReaction:
            if not inputReaction[i]['kegg']=='C00001' and inputReaction[i]['conc']:
                dG_correction += inputReaction[i]['stochio']*np.log(inputReaction[i]['conc'])
        return dG0_r_prime+(dG_correction*self.RT), dG0_uncertainty




    #compound is a dictionnary with a KEGG within compound_dict and 'stochio'
    def compound_to_matrix(self, compound):
        x = np.array(np.zeros((np.array(self.cc_preprocess['C1']).shape[0], 1)), ndmin=2)
        g = np.array(np.zeros((np.array(self.cc_preprocess['C3']).shape[0], 1)), ndmin=2)
        # x is the stoichiometric vector of the reaction, only for the
        # compounds that appeared in the original training set for CC
        compound_index = self.compound_dict[compound['kegg']]['compound_index']
        if compound_index is None:
            logging.error('could not find index for '+str(cmp_name))
            return x, g
        x[compound_index, 0] = 1
        x = compound['stochio']*x
        # g is the group incidence vector of all the other compounds
        gv = self.compound_dict[compound['kegg']]['group_vector']
        if gv is None:
            logging.error('could not find group vector for '+str(cmp_name))
        for g_ind, g_count in gv:
            g[g_ind, 0] += g_count
        g = compound['stochio']*g
        return x, g
    '''


    '''
    ##might be able to remove all these functions and replace them with single components as we only need to perform this calcuation one compound at a time
    def reaction_dG0_prime(self, inputReaction):
        """
            Calculate the standard dG'0 of reaction.
            Arguments:
        """
        dG0_r_prime = 0
        for cmp_name in inputReaction:
            if inputReaction[cmp_name]['kegg']:
                dG0_f_prime = compound_dG0(inputReaction[cmp_name]['kegg'])
                if dG0_f_prime is None:
                    return None
                dG0_r_prime += inputReaction[cmp_name]['stochio']*dG0_f_prime
        return dG0_r_prime

    def reaction_to_matrix(self, inputReaction, compound_dict):
        Nc = np.array(self.cc_preprocess['C1']).shape[0]
        Ng = np.array(self.cc_preprocess['C3']).shape[0]
        x = np.array(np.zeros((Nc, 1)), ndmin=2)
        g = np.array(np.zeros((Ng, 1)), ndmin=2)
        for cmp_name in inputReaction:
            # x is the stoichiometric vector of the reaction, only for the
            # compounds that appeared in the original training set for CC
            x_c = np.array(np.zeros((Nc, 1)))
            compound_index = compound_dict[inputReaction[cmp_name]['kegg']]['compound_index']
            if compound_index is None:
                logging.error('could not find index for '+str(cmp_name))
            x_c[compound_index, 0] = 1
            x += inputReaction[cmp_name]['stochio']*x_c
            # g is the group incidence vector of all the other compounds
            g_c = np.array(np.zeros((Ng, 1)))
            gv = compound_dict[inputReaction[cmp_name]['kegg']]['group_vector']
            if gv is None:
                logging.error('could not find group vector for '+str(cmp_name))
            for g_ind, g_count in gv:
                g_c[g_ind, 0] += g_count
            g += inputReaction[cmp_name]['stochio']*g_c
        return x, g
    '''














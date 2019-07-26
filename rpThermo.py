import numpy as np
import logging
import pybel
import json
from scipy.special import logsumexp
import copy
import pickle
import os
import libsbml

#local package
#import component_contribution
from . import component_contribution

class rpThermo:
    """Combination of equilibrator and group_contribution analysis to calculate the thermodymaics of the individual
    metabolic pathways output from RP2paths and thus RetroPath2.0
    """
    def __init__(self, 
            pH=7.0, 
            pMg=14.0, 
            ionic_strength=0.1, 
            temperature=298.15):
        #self.compound_dict = compound_dict
        #TODO: put this in its own function
        #### set the physiological parameters
        self.pH = pH
        self.pMg = pMg
        self.ionic_strength = ionic_strength
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
        self.debye_huckel = self.DH_alpha*self.ionic_strength**(0.5)/(1.0+self.DH_beta*self.ionic_strength**(0.5))
        #temporarely save the calculated dG's
        self.calculated_dG = {}
        if not self._loadCache():
            raise ValueError
        

    ################ PRIVATE FUNCTIONS ###################

    
    ##
    #
    #
    def _loadCache(self):
        dirname = os.path.dirname(os.path.abspath( __file__ ))
        try:
            self.cc_preprocess = np.load(dirname+'/cache/cc_preprocess.npz')
        except FileNotFoundError as e:
            logging.error(e)
            return False
        try:
            self.kegg_dG = pickle.load(open(dirname+'/cache/kegg_dG.pickle', 'rb'))
        except FileNotFoundError as e:
            logging.error(e)
            return False
        return True


    ## Select a thermodynamics score from the dataset
    #
    # Given that there can be multiple precalculated dG for a given molecule (MNXM)
    # this function gives the priority to a particular order for a given MNXM ID
    # alberty>smallest(KEGG_ID)>other(KEGG_ID)
    #
    def _select_dG(self, cid):
        if cid in self.kegg_dG:
            if 'component_contribution' in self.kegg_dG[cid] and self.kegg_dG[cid]['component_contribution']:
                ####### select the smallest one
                '''
                if len(self.kegg_dG[cid]['component_contribution'])==1:
                    return self.kegg_dG[cid]['component_contribution'][0]['compound_index'], self.kegg_dG[cid]['component_contribution'][0]['group_vector'], self.kegg_dG[cid]['component_contribution'][0]['pmap']['species']
                else:
                    #select the lowest CID and make sure that there
                    toRet = {'CID': 'C99999'}
                    for cmp_dict in self.kegg_dG[cid]['component_contribution']:
                        if int(cmp_dict['CID'][1:])<int(toRet['CID'][1:]):
                            toRet = cmp_dict
                    return toRet['compound_index'], toRet['group_vector'], toRet['pmap']['species']
                '''
                ###### select the one with the most information
                #return the smallest one that has all the information required
                for cmp_d in sorted(self.kegg_dG[cid]['component_contribution'], key=lambda k: int(k['CID'][1:])):
                    if 'compound_index' in cmp_d and 'group_vector' in cmp_d and 'pmap' in cmp_d and 'species' in cmp_d['pmap']:
                        return cmp_d['compound_index'], cmp_d['group_vector'], cmp_d['pmap']['species']
                #if cannot find return the smallest one and return the info that you can
                compound_index = None
                group_vector = None
                pmap_species = None
                try:
                    compound_index = self.kegg_dG[cid]['component_contribution'][0]['compound_index'] 
                except KeyError:
                    pass
                try:
                    group_vector = self.kegg_dG[cid]['component_contribution'][0]['group_vector']
                except KeyError:
                    pass
                try:
                    pmap_species = self.kegg_dG[cid]['component_contribution'][0]['pmap']['species']
                except KeyError:
                    logging.warning('Component_contribution cannot find any species')
                    raise KeyError
                #if compound_index==None and group_vector==None and pmap_species==None:
                    #exit the component_conribution condition if all are none
                #    continue
                return compound_index, group_vector, pmap_species
            #if you cannot find the component in the component_contribution then select the alberty dataset
            elif 'alberty' in self.kegg_dG[cid] and self.kegg_dG[cid]['alberty']:
                if len(self.kegg_dG[cid]['alberty'])==1:
                    return None, None, self.kegg_dG[cid]['alberty'][0]['species']
                else:
                    logging.warning('Alberty and component_contribution species are empty')
                    raise KeyError
            else:
                logging.warning('There are no valid dictionnary of precalculated dG for '+str(cid))
                raise KeyError
        else:
            logging.warning('There are no '+str(cid)+' in self.kegg_dG')
            raise KeyError


    ####################################################
    ################# ON THE FLY #######################
    ####################################################


    #must make sure that there is no KEGG id and that there is a SMILES description
    #TODO export the matrices
    #TODO add the option of accepting InChI strings as well as SMILES


    ## Calculate the formation energy of a compound
    #
    # Script adapted from the Noor et al.'s component contribution project. Seperates a compound into groups
    # Decompose a SMILES string
    # Warning -- the dimensions of x and g are not the same as the compound_to_matrix function
    # calculate pKas of the target and intermediates using cxcalc
    # and calculates the its formation energy from the individual contribution of each sub-compound.
    #
    def scrt_dfG_prime_o(self, srct_type, srct_string, stoichio):
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
        X = np.zeros((self.cc_preprocess['C1'].shape[0], 1))
        G = np.zeros((self.cc_preprocess['C3'].shape[0], 1))
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
        elif 0<self.ionic_strength:
            ddG = sum(p_kas[0:major_microspecies]) * self.R * self.temperature * np.log(10)
        else:
            ddG = -sum(p_kas[major_microspecies:0]) * self.R * self.temperature * np.log(10)
        ddG0_forward = trans+ddG
        dG0 = dG0_cc+ddG0_forward
        dG0_prime = stoichio*dG0
        if type(dG0_prime)==np.ndarray:
            dG0_prime = dG0_prime[0].astype(float)
        return dG0_prime, X, G, {'atom_bag': atom_bag, 'p_kas': p_kas, 'major_ms': major_microspecies, 'number_of_protons': number_of_protons, 'charges': charges}


    ##########################################################
    ################### PRECALCULATED ########################
    ##########################################################


    ## Calculate the formation energy taking into consideration pH and ionic strength of environment
    #from the formation energy of the compound calculate its prime
    #ie. taking into account the pH and the ionic strength of the environment
    def dG_dGprime(self, nMg, z, nH, dG0_f):
        sqrt_I = np.sqrt(self.ionic_strength)
        dG0_prime = dG0_f
        # add the potential related to the pH
        if nH > 0:
            dG0_prime += nH*self.RTlog10*self.pH
        # add the potential related to the ionic strength
        dG0_prime -= self.debye_hueckle_a*(z**2-nH)*sqrt_I/(1+self.debye_hueckle_b*sqrt_I)
        # add the potential related to the Mg ions
        if nMg > 0:
            dG0_prime += nMg*(self.RTlog10*nMg-self.mg_formation_energy)
        #print('nH = %d, z = %d, dG0 = %.1f --> dG0\' = %.1f' %
        #              (nH, z, dG0_f, dG0_prime))
        return dG0_prime


    ## Calculate the formation energy of a compoud at 1mM
    #
    #
    def cmp_dfG_prime_o(self, kegg_cid, stoichio):
        ########### dG0_f ###########
        #WARNING: we are ignoring gas phase even if it would be valid in some cases (O2 for example)
        physioParameter = None #determine the phase for the concentration adjustement 
        try:
            compound_index, group_vector, species_list = self._select_dG(kegg_cid)
            scaled_transforms = [-self.dG_dGprime(i['nMg'], 
                                           i['z'], 
                                           i['nH'], 
                                           i['dG0_f'])/self.RT for i in species_list if not i['phase']=='gas']
            #note that this should change if a user wants to input his own values....
            #WARNING: we are taking the first one only
            if species_list[0]['phase']=='aqueous':
                physioParameter = 1e-3
            else:
                physioParameter = 1.0
        except KeyError:
            logging.warning('Cannot find precalculated cmp_dfG_prime_o')
            raise KeyError
        total = scaled_transforms[0]
        for i in range(1, len(scaled_transforms)):
            total = np.logaddexp(total, scaled_transforms[i])
        ############ uncertainty #############
        X = np.zeros((self.cc_preprocess['C1'].shape[0], 1))
        G = np.zeros((self.cc_preprocess['C3'].shape[0], 1))
        if compound_index is not None:
            X[compound_index, 0] = stoichio
        elif group_vector is not None:
            for g_ind, g_count in group_vector:
                G[g_ind, 0] += stoichio*g_count
        else:
            logging.warning('Cannot generate uncerstainty for '+str(kegg_cid))
        dG0_prime = -self.RT*total
        return dG0_prime, X, G, physioParameter


    ##########################################################
    ################## BOTH ##################################
    ##########################################################
    

    ## Calculate the uncertainty with a given Gibbs free enery
    #
    #
    def dG0_uncertainty(self, X, G):
        return 1.92*float(np.sqrt(X.T @ self.cc_preprocess['C1'] @X +
                             X.T @ self.cc_preprocess['C2'] @G +
                             G.T @ self.cc_preprocess['C3'] @G ))

   
    ## takes a list of SBase libsbml objects and extracts the stochiometry from it
    #
    #TODO: extract the concentration (if defined) instead of using the milliMolar concentration adjustment
    def concentrationCorrection(self, stoichio, conc=None):
        #if no concentrations are specified then we assume millimolar adjustment
        #return self.R*self.temperature*(stoichio*np.log(conc))
        if conc==None:
            conc = [1e-3]*len(stochio)
        return self.R*self.temperature*sum([sto*np.log(co) for sto, co in zip(stoichio, conc)])


    ###########################################################
    ################# RPmodel update ##########################
    ###########################################################


    ## Calculate a  species dG0_prime_o and its uncertainty
    #
    #
    def species_dfG_prime_o(self, species_annot, stoichio):
        #check to see if there are mutliple, non-deprecated, MNX ids in annotations
        X = None
        G = None
        dfG_prime_o = None
        cid = None
        physioParameter = None #this paraemter determines the concentration of the copound for the adjustemet, (dG_prime_m). It assumes physiological conditions ie 1e-3 for aquaeus and 1 for gas, solid etc.... TODO: Next step is to have the user input his own
        #Try to find your species in the already calculated species
        smiles = species_annot.getChild('RDF').getChild('Ibisba').getChild('ibisba').getChild('smiles').getChild(0).toXMLString()
        inchi = species_annot.getChild('RDF').getChild('Ibisba').getChild('ibisba').getChild('inchi').getChild(0).toXMLString()
        ################### use already calculated ##########################
        if inchi in self.calculated_dG:
            X = self.calculated_dG[inchi]['X']
            G = self.calculated_dG[inchi]['G']
            dfG_prime_o = self.calculated_dG[inchi]['dfG_prime_o']
        elif smiles in self.calculated_dG:
            X = self.calculated_dG[smiles]['X']
            G = self.calculated_dG[smiles]['G']
            dfG_prime_o = self.calculated_dG[smiles]['dfG_prime_o']
        else:
            ################## use the KEGG id for precalculated ######################
            #if not try to find it in cc_preprocess
            cids = []
            bag = species_annot.getChild('RDF').getChild('Description').getChild('is').getChild('Bag')
            for i in range(bag.getNumChildren()):
                str_annot = bag.getChild(i).getAttrValue(0)
                if str_annot.split('/')[-2]=='kegg.compound':
                    cids.append(str_annot.split('/')[-1])
            cid = [i for i in cids if i[0]=='C'] #some of the KEGG compounds are drugs and start with a D
            if len(cid)==1:
                cid = cid[0]
            #TODO: sort them to have the lowest KEGG CID used
            elif len(cid)>1:
                cid = sorted(cid)[0]
                logging.warning('There are more than one results, using the first one '+str(sorted(cid)))
            else:
                cid = None
            if cid:
                #BUG or overlooked? overwrite H+ and H2O dfG_prime_m
                if not cid=='C00080': #ignore H+ for the moment
                    try:
                        dfG_prime_o, X, G, physioParameter = self.cmp_dfG_prime_o(
                                cid, 
                                stoichio)
                    except KeyError:
                        #TODO: seperate this part in a function instead of repeating the code
                        #calculate using the structure
                        #try for smiles
                        if inchi:
                            try:
                                dfG_prime_o, X, G, additional_info = self.scrt_dfG_prime_o(
                                        'inchi', 
                                        inchi, 
                                        stoichio)
                                self.calculated_dG[smiles] = {}
                                self.calculated_dG[smiles]['dfG_prime_o'] = dfG_prime_o
                                self.calculated_dG[smiles]['X'] = X
                                self.calculated_dG[smiles]['G'] = G
                            except (KeyError, LookupError):
                                logging.warning('Cannot use InChI to calculate the thermodynamics')
                                pass
                        #try for inchi
                        if smiles and dfG_prime_o==None:
                            try:
                                dfG_prime_o, X, G, additional_info = self.scrt_dfG_prime_o(
                                        'smiles',
                                        smiles,
                                        stoichio)
                                self.calculated_dG[smiles] = {}
                                self.calculated_dG[smiles]['dfG_prime_o'] = dfG_prime_o
                                self.calculated_dG[smiles]['X'] = X
                                self.calculated_dG[smiles]['G'] = G
                            except (KeyError, LookupError):
                                logging.warning('Cannot use SMILES to calculate the thermodynamics')
                                raise KeyError
            ############## use the structure #########################
            else:
                if cid:
                    logging.warning('Database does not contain thermodynamics for KEGG: '+str(cid))
                #at last, if all fails use its structure to calculate dfG_prime_o
                #try for smiles
                if inchi:
                    try:
                        dfG_prime_o, X, G, additional_info = self.scrt_dfG_prime_o(
                                'inchi', 
                                inchi, 
                                stoichio)
                        self.calculated_dG[smiles] = {}
                        self.calculated_dG[smiles]['dfG_prime_o'] = dfG_prime_o
                        self.calculated_dG[smiles]['X'] = X
                        self.calculated_dG[smiles]['G'] = G
                    except (KeyError, LookupError):
                        logging.warning('Cannot use InChI to calculate the thermodynamics')
                        pass
                #try for inchi
                if smiles and dfG_prime_o==None:
                    try:
                        dfG_prime_o, X, G, additional_info = self.scrt_dfG_prime_o(
                                'smiles',
                                smiles,
                                stoichio)
                        self.calculated_dG[smiles] = {}
                        self.calculated_dG[smiles]['dfG_prime_o'] = dfG_prime_o
                        self.calculated_dG[smiles]['X'] = X
                        self.calculated_dG[smiles]['G'] = G
                    except (KeyError, LookupError):
                        logging.warning('Cannot use SMILES to calculate the thermodynamics')
                        raise KeyError
        #update the species dfG information
        dfG_prime_m = None
        if physioParameter==None:
            physioParameter = 1e-3
        if not dfG_prime_o==None:
            write_dfG_prime_o = dfG_prime_o
            dfG_prime_m = dfG_prime_o+self.concentrationCorrection([abs(stoichio)], [physioParameter]) 
            write_dfG_prime_m = dfG_prime_m
            write_uncertainty = self.dG0_uncertainty(X, G)
        elif cid=='C00080':
            write_dfG_prime_o = 0.0
            write_dfG_prime_m = -17.1
            write_uncertainty = 5.8
            ### this is wrong.... need to define it better
            X = np.zeros((self.cc_preprocess['C1'].shape[0], 1))
            G = np.zeros((self.cc_preprocess['C3'].shape[0], 1)) 
        else:
            write_dfG_prime_o = 0.0
            write_dfG_prime_m = 0.0
            write_uncertainty = 0.0
            X = np.zeros((self.cc_preprocess['C1'].shape[0], 1))
            G = np.zeros((self.cc_preprocess['C3'].shape[0], 1)) 
        ibisba_annot = species_annot.getChild('RDF').getChild('Ibisba').getChild('ibisba')
        tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<ibisba:ibisba xmlns:ibisba="http://ibisba.eu"> <ibisba:dfG_prime_o units="kj_per_mol" value="'+str(write_dfG_prime_o)+'" /> </ibisba:ibisba>')
        ibisba_annot.addChild(tmpAnnot.getChild('dfG_prime_o'))
        tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<ibisba:ibisba xmlns:ibisba="http://ibisba.eu"> <ibisba:dfG_prime_m units="kj_per_mol" value="'+str(write_dfG_prime_m)+'" /> </ibisba:ibisba>')
        ibisba_annot.addChild(tmpAnnot.getChild('dfG_prime_m'))
        tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<ibisba:ibisba xmlns:ibisba="http://ibisba.eu"> <ibisba:dfG_uncert units="kj_per_mol" value="'+str(write_uncertainty)+'" /> </ibisba:ibisba>')
        ibisba_annot.addChild(tmpAnnot.getChild('dfG_uncert'))
        return dfG_prime_o, X, G, physioParameter #we call physioParameter concentration later 


    ## Calculate the pathway and reactions dG0_prime_o and its uncertainty
    #
    # WARNING: we skip the components that we fail to calculate the thermodynamics and 
    # and as a consequence the pathway thermo score ignore that step --> need to at least report it
    def pathway_drG_prime_m(self, rpsbml, pathId='rp_pathway'):
        #calculate the number of species that are involded in the pathway 
        #for each species calculate the ddG_prime_m and its uncertainty
        #test to see if there are mutliple MNX and remove all deprecated ones
        groups = rpsbml.model.getPlugin('groups')
        rp_pathway = groups.getGroup(pathId)
        #get the number of species, reactions and pathways ###
        pathway_dfG_prime_o = 0.0
        X_path = np.zeros((self.cc_preprocess['C1'].shape[0], 1))
        G_path = np.zeros((self.cc_preprocess['C3'].shape[0], 1))
        pathway_stoichio = []
        pathway_concentration = []
        #path
        already_calculated = {} 
        for member in rp_pathway.getListOfMembers():
            reaction_dfG_prime_o = 0.0
            X_reaction = np.zeros((self.cc_preprocess['C1'].shape[0], 1))
            G_reaction = np.zeros((self.cc_preprocess['C3'].shape[0], 1))
            reaction_stoichio = []
            reaction_concentration = []
            reaction = rpsbml.model.getReaction(member.getIdRef())
            #react
            for pro in reaction.getListOfProducts():
                #rpsbml.model.getSpecies(pro.species).getAnnotation()
                #print(rpsbml.model.getSpecies(pro.species).getAnnotation().toXMLString())
                try:
                    if not pro.species in already_calculated:
                        dfG_prime_o, X, G, concentration = self.species_dfG_prime_o(
                                rpsbml.model.getSpecies(pro.species).getAnnotation(), 
                                float(pro.stoichiometry))
                        already_calculated[pro.species] = {}
                        already_calculated[pro.species]['dfG_prime_o'] = dfG_prime_o
                        already_calculated[pro.species]['X'] = X
                        already_calculated[pro.species]['G'] = G
                        already_calculated[pro.species]['concentration'] = concentration
                    else:
                        dfG_prime_o = already_calculated[pro.species]['dfG_prime_o']
                        X = already_calculated[pro.species]['X']
                        G = already_calculated[pro.species]['G']
                        concentration = already_calculated[pro.species]['concentration']
                except (KeyError, LookupError):
                    logging.warning('Failed to calculate the thermodynamics for '+str(pro.species))
                    continue
                reaction_stoichio.append(float(pro.stoichiometry))
                pathway_stoichio.append(float(pro.stoichiometry))
                reaction_concentration.append(float(concentration))
                pathway_concentration.append(float(concentration))
                if not dfG_prime_o==None:
                    reaction_dfG_prime_o += dfG_prime_o*pro.stoichiometry
                    pathway_dfG_prime_o += dfG_prime_o*pro.stoichiometry
                    X_reaction += X
                    G_reaction += G
                    X_path += X
                    G_path += G
            for rea in reaction.getListOfReactants():
                #print(rpsbml.model.getSpecies(rea.species).getAnnotation().toXMLString())
                try:
                    if not rea.species in already_calculated:
                        dfG_prime_o, X, G, concentration = self.species_dfG_prime_o(
                                rpsbml.model.getSpecies(rea.species).getAnnotation(),
                                -float(rea.stoichiometry))
                        already_calculated[rea.species] = {}
                        already_calculated[rea.species]['dfG_prime_o'] = dfG_prime_o
                        already_calculated[rea.species]['X'] = X
                        already_calculated[rea.species]['G'] = G
                        already_calculated[rea.species]['concentration'] = concentration
                    else:
                        dfG_prime_o = already_calculated[rea.species]['dfG_prime_o']
                        X = already_calculated[rea.species]['X']
                        G = already_calculated[rea.species]['G']
                        concentration = already_calculated[rea.species]['concentration']
                except (KeyError, LookupError):
                    logging.warning('Failed to calculate the thermodynamics for '+str(pro.species))
                    continue
                reaction_stoichio.append(-float(rea.stoichiometry))
                pathway_stoichio.append(-float(rea.stoichiometry))
                reaction_concentration.append(float(concentration))
                pathway_concentration.append(float(concentration))
                if not dfG_prime_o==None:
                    reaction_dfG_prime_o += dfG_prime_o*(-rea.stoichiometry)
                    pathway_dfG_prime_o += dfG_prime_o*(-rea.stoichiometry)
                    X_reaction += X
                    G_reaction += G
                    X_path += X
                    G_path += G
            #add the reaction thermo to the sbml
            #ibisba_annot = member.getIdRef().getAnnotation().getChild('RDF').getChild('Ibisba').getChild('ibisba')
            reac = rpsbml.model.getReaction(member.getIdRef())
            reac_annot = reac.getAnnotation()
            ibisba_annot = reac_annot.getChild('RDF').getChild('Ibisba').getChild('ibisba')
            tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<ibisba:ibisba xmlns:ibisba="http://ibisba.eu"> <ibisba:dfG_prime_o units="kj_per_mol" value="'+str(reaction_dfG_prime_o)+'" /> </ibisba:ibisba>')
            ibisba_annot.addChild(tmpAnnot.getChild('dfG_prime_o'))
            tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<ibisba:ibisba xmlns:ibisba="http://ibisba.eu"> <ibisba:dfG_prime_m units="kj_per_mol" value="'+str(reaction_dfG_prime_o+self.concentrationCorrection(reaction_stoichio, reaction_concentration))+'" /> </ibisba:ibisba>')
            ibisba_annot.addChild(tmpAnnot.getChild('dfG_prime_m'))
            tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<ibisba:ibisba xmlns:ibisba="http://ibisba.eu"> <ibisba:dfG_uncert units="kj_per_mol" value="'+str(self.dG0_uncertainty(X_reaction, G_reaction))+'" /> </ibisba:ibisba>')
            ibisba_annot.addChild(tmpAnnot.getChild('dfG_uncert'))
            '''
            rpsbml.model.getReaction(member.getIdRef()).getAnnotation().getChild('RDF').getChild('Ibisba').getChild('ibisba').getChild('dG_prime_o').addAttr('value', str(reaction_dfG_prime_o))
            rpsbml.model.getReaction(member.getIdRef()).getAnnotation().getChild('RDF').getChild('Ibisba').getChild('ibisba').getChild('dG_prime_o').addAttr('value', str(reaction_dfG_prime_o))
            #rpsbml.model.getReaction(member.getIdRef()).getAnnotation().getChild('RDF').getChild('Ibisba').getChild('ibisba').getChild('dG_prime_m').addAttr('value', str(reaction_dfG_prime_m))
            rpsbml.model.getReaction(member.getIdRef()).getAnnotation().getChild('RDF').getChild('Ibisba').getChild('ibisba').getChild('dG_prime_m').addAttr('value', str(reaction_dfG_prime_o+self.concentrationCorrection(reaction_stoichio, reaction_concentration)))
            rpsbml.model.getReaction(member.getIdRef()).getAnnotation().getChild('RDF').getChild('Ibisba').getChild('ibisba').getChild('dG_uncert').addAttr('value', str(self.dG0_uncertainty(X_reaction, G_reaction)))
            '''
        #add the pathway thermo to the sbml
        ibisba_annot = rp_pathway.getAnnotation().getChild('RDF').getChild('Ibisba').getChild('ibisba')
        tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<ibisba:ibisba xmlns:ibisba="http://ibisba.eu"> <ibisba:dfG_prime_o units="kj_per_mol" value="'+str(pathway_dfG_prime_o)+'" /> </ibisba:ibisba>')
        ibisba_annot.addChild(tmpAnnot.getChild('dfG_prime_o'))
        tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<ibisba:ibisba xmlns:ibisba="http://ibisba.eu"> <ibisba:dfG_prime_m units="kj_per_mol" value="'+str(pathway_dfG_prime_o+self.concentrationCorrection(pathway_stoichio, pathway_concentration))+'" /> </ibisba:ibisba>')
        ibisba_annot.addChild(tmpAnnot.getChild('dfG_prime_m'))
        tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<ibisba:ibisba xmlns:ibisba="http://ibisba.eu"> <ibisba:dfG_uncert units="kj_per_mol" value="'+str(self.dG0_uncertainty(X_path, G_path))+'" /> </ibisba:ibisba>')
        ibisba_annot.addChild(tmpAnnot.getChild('dfG_uncert'))
        '''
        rp_pathway.getAnnotation().getChild('RDF').getChild('Ibisba').getChild('ibisba').getChild('dG_prime_o').addAttr('value', str(pathway_dfG_prime_o))
        #rp_pathway.getAnnotation().getChild('RDF').getChild('Ibisba').getChild('ibisba').getChild('dG_prime_m').addAttr('value', str(pathway_dfG_prime_m))
        rp_pathway.getAnnotation().getChild('RDF').getChild('Ibisba').getChild('ibisba').getChild('dG_prime_m').addAttr('value', str(pathway_dfG_prime_o+self.concentrationCorrection(pathway_stoichio, pathway_concentration)))
        rp_pathway.getAnnotation().getChild('RDF').getChild('Ibisba').getChild('ibisba').getChild('dG_uncert').addAttr('value', str(self.dG0_uncertainty(X_path, G_path)))
        '''


    #TODO: implement to return if the reaction is balanced and the reversibility index 

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

if __name__ == "__main__":
    #READ THE TAR.XZ FILE
    rpsbml_paths = {}
    tar = tarfile.open('tests/testTHERMOin.tar.xz') #TODO: create this
    rpsbml_paths = {}
    for member in tar.getmembers():
        rpsbml_paths[member.name] = rpFBA.rpSBML(member.name,libsbml.readSBMLFromString(tar.extractfile(member).read().decode("utf-8")))
    ###
    rpthermo = rpThermo()
    for rpsbml_name in rpsbml_paths:
        rpthermo.pathway_drG_prime_m(rpsbml_paths[rpsbml_name])
    #WRITE THE TAR.XZ FILE
    with tarfile.open('testFBAout.tar.xz', 'w:xz') as tf:
        for rpsbml_name in rpsbml_paths:
            data = libsbml.writeSBMLToString(rpsbml_paths[rpsbml_name].document).encode('utf-8')
            fiOut = BytesIO(data)
            info = tarfile.TarInfo(rpsbml_name)
            info.size = len(data)
            tf.addfile(tarinfo=info, fileobj=fiOut)

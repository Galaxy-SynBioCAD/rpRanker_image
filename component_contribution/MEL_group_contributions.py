#Calculate the formation energy using the group contributions only (basically not using the KEGG identifiers)

reac_smiles = []
params = dict(np.load('/home/mdulac/workspace/component-contribution/cache/component_contribution.npz'))
params = dict(np.load('~/workspace/equilibrator-api/equilibrator_api/data/cc_preprocess.npz'))
G1 = params['preprocess_G1']
G2 = params['preprocess_G2']
G3 = params['preprocess_G3']
groups_data = component_contribution.inchi2gv.init_groups_data()
group_names = groups_data.GetGroupNames()
decomposer = component_contribution.inchi2gv.InChIDecomposer(groups_data)

##################### dG0_r ###########################
g = np.zeroes(len(group_names)) #array for the same length of the groups that we will decompose the reaction with
for smiles in reac_smiles:
    #try: #try without this to see what are the returned errors to catch
    g += decomposer.smiles_to_groupvec(smiles) 
    #except inchi2gv.GroupDecompositionError as exception:
    #    print('Error decomposing the SMILES: '+str(smiles))
    #    raise exception
weights = (g @ G3[0:len(group_names), :]).round(5)
orders = sorted(range(weights.shape[1]), key=lambda j:abs(weights[0, j]), reverse=True)
# dG0_cc = (x*G1 + x*G2 + g*G3)*b

##################### dG0_r_prime #####################
#dG0_prime = dG0_cc + reaction.get_transform_ddG0(pH, I, T) #requires the compound from cache


# In the program -- looks for a compound and if it does not find it then queries cxcalc -- we use cxcalc all the time
cpd = Compound.get(compound_id, compute_pkas)

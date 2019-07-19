#list of functions that are usefull but not directly related to the rpFBA
#import cobra
import csv
#import libsbml
import rpSBML
import pickle
import gzip
import os

## Generate the sink from a given model and the 
#
# NOTE: this only works for MNX models, since we are parsing the id
# TODO: change this to read the annotations and extract the MNX id's
#
def genSink(sbml_model, file_out, compartment_id='MNXC3'):
    ### open the cache ###
    dirname = os.path.dirname(os.path.abspath( __file__ ))
    smiles_inchi = pickle.load(gzip.open(dirname+'/cache/smiles_inchi.pickle.gz', 'rb'))
    rpsbml = rpSBML.rpSBML('tmp')
    rpsbml.readModel(sbml_model)
    cytoplasm_species = []
    for i in rpsbml.model.getListOfSpecies():
        if i.getCompartment()==compartment_id:
            cytoplasm_species.append(i)
    with open(file_out, mode='w') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
        writer.writerow(['Name','InChI'])
        for i in cytoplasm_species:
            res = rpsbml.readAnnotation(i.getAnnotation())
            #extract the MNX id's
            try:
                mnx = res['metanetx.chemical'][0]
            except KeyError:
                continue
            #mnx = i.getId().split('__')[0]
            try:
                inchi = meta_inchi[mnx]
            except KeyError:
                inchi = None
            if mnx and inchi:
                writer.writerow([mnx,inchi])

'''
def genSink(path_cobraModel, path_chem_prop, file_out_name=None, compartment=None):
    """Function to extract from an SBML MNX model the available compounds
    and generating the sink for RetroPath2.0. There is the option to
    generate the sink from a particular model compartment
    """
    #open the model file and extract all the metabolites; make sure they are unique
    model = cobra.io.read_sbml_model(path_cobraModel)
    if compartment:
        meta = set([i.id.split('__')[0] for i in model.metabolites if i.id.split('__')[2]==compartment])
    else:
        meta = set([i.id.split('__')[0] for i in model.metabolites])
    #open the MNX list of all chemical species
    meta_inchi = {}
    with open(path_chem_prop) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if not row[0][0]=='#' and row[0][:3]=='MNX' and not row[5]=='':
                meta_inchi[row[0]] = row[5]
    #write the results to a new file
    if not file_out_name:
        file_out_name = 'sink.csv'
    with open(file_out_name, mode='w') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
        writer.writerow(['Name','InChI'])
        for m in meta:
            try:
                writer.writerow([m, meta_inchi[m]])
            except KeyError:
                print('Cannot find '+str(m)+' in '+str(path_chem_prop))
'''

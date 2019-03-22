#list of functions that are usefull but not directly related to the rpFBA
import cobra
import csv

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

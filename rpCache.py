import csv
import logging
import gzip
import os
import json
import pickle
import sqlite3
import re
import itertools
from ast import literal_eval
import gzip
from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi, MolToInchiKey, AddHs
from shutil import copyfile

## @package Cache
#
# Documentation for the cache generation of rpFBA


## \brief Class to generate the cache
#
# Contains all the functions that parse different files, used to calculate the thermodynamics and the FBA of the 
#the other steps. These should be called only when the files have changes
class rpCache:
    ## Cache constructor
    # 
    # @param self The object pointer
    # @param inputPath The path to the folder that contains all the input/output files required
    def __init__(self):
        #given by Thomas
        self.convertMNXM = {'MNXM162231': 'MNXM6',
                'MNXM84': 'MNXM15',
                'MNXM96410': 'MNXM14',
                'MNXM114062': 'MNXM3',
                'MNXM145523': 'MNXM57',
                'MNXM57425': 'MNXM9',
                'MNXM137': 'MNXM588022'}
        #personally looked at the KEGG to MNXM conversion for the thermodynamics 
        self.deprecatedMNXM_mnxm = None

    #######################################################
    ################### PRIVATE FUNCTION ##################
    #######################################################



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
    def _convert_depiction(self, idepic, itype='smiles', otype={'inchikey'}):
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



    ########################################################
    ####################### PUBLIC FUNCTIONS ###############
    ######################################################## 


    #[TODO] merge the two functions
    ## Function to parse the chem_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of old to new MetanetX identifiers
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return a The dictionnary of identifiers  
    #TODO: save the self.deprecatedMNXM_mnxm to be used in case there rp_paths uses an old version of MNX
    def deprecatedMNXM(self, chem_xref_path):
        a = {}
        with open(chem_xref_path) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                mnx = row[0].split(':')
                if not row[0][0]=='#':
                    if mnx[0]=='deprecated':
                        #a[row[1]] = mnx[1]
                        a[mnx[1]] = row[1]
            a.update(self.convertMNXM)
            a['MNXM01'] = 'MNXM1'
        return a

    ## Function to parse the chem_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of old to new MetanetX identifiers
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return a The dictionnary of identifiers  
    #TODO: save the self.deprecatedMNXM_mnxm to be used in case there rp_paths uses an old version of MNX
    def chemXref(self, chem_xref_path):
        chemXref = {}
        with open(chem_xref_path) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnx = row[1]
                    if len(row[0].split(':'))==1:
                        dbName = 'mnx'
                        dbId = row[0]
                    else:
                        dbName = row[0].split(':')[0]
                        dbId = ''.join(row[0].split(':')[1:])
                        if dbName=='deprecated':
                            dbName = 'mnx'
                    ### MNXM ###
                    if not mnx in chemXref:
                        chemXref[mnx] = {}
                    if not dbName in chemXref[mnx]:
                        chemXref[mnx][dbName] = []
                    if not dbId in chemXref[mnx][dbName]:
                        chemXref[mnx][dbName].append(dbId)
                    ### DB ###
                    if not dbName in chemXref:
                        chemXref[dbName] = {}
                    if not dbId in chemXref[dbName]:
                        chemXref[dbName][dbId] = mnx
        return chemXref


    ## Function to parse the reacXref.tsv file of MetanetX
    #
    #  Generate a dictionnary of old to new MetanetX identifiers
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return a The dictionnary of identifiers  
    #TODO: save the self.deprecatedMNXM_mnxm to be used in case there rp_paths uses an old version of MNX
    def reacXref(self, reacXref_path, reac_prop_path):
        reacXref = {}
        with open(reacXref_path) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#' and len(row[0].split(':'))==2:
                    mnx = row[1]
                    if len(row[0].split(':'))==1:
                        dbName = 'mnx'
                        dbId = row[0]
                    else:
                        dbName = row[0].split(':')[0]
                        dbId = ''.join(row[0].split(':')[1:])
                        if dbName=='deprecated':
                            dbName = 'mnx'
                    if not mnx in reacXref:
                        reacXref[mnx] = {}
                    if not dbName in reacXref[mnx]:
                        reacXref[mnx][dbName] = []
                    if not dbId in reacXref[mnx][dbName]:
                        reacXref[mnx][dbName].append(dbId)
        #use this to retreive the EC number for the reactions
        with open(reac_prop_path) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#' and not row[4]=='':
                    mnx = row[0]
                    if not mnx in reacXref:
                        reacXref[mnx] = {}
                    dbName = 'ec'
                    if not dbName in reacXref[mnx]:
                        reacXref[mnx][dbName] = []
                    for ec in row[4].split(';'):
                        if not ec in reacXref[mnx][dbName]:
                            reacXref[mnx][dbName].append(ec)
        return reacXref



    ## Function to parse the compXref.tsv file of MetanetX
    #
    #  Generate a dictionnary of old to new MetanetX identifiers
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return a The dictionnary of identifiers  
    #TODO: save the self.deprecatedMNXM_mnxm to be used in case there rp_paths uses an old version of MNX
    def compXref(self, compXref_path, comp_prop_path):
        mnxc_name = {}
        with open(comp_prop_path) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    if not row[0] in mnxc_name:
                        mnxc_name[row[0]] = row[1]
        #Mel --> not a fan of the hardcoding, if the file changes then one would need to add new entries
        #possCID = ['mnxc', 'bigg', 'cco', 'go', 'seed', 'name']
        #pubDB_name_xref = {}
        name_pubDB_xref = {}
        try:
            with open(compXref_path) as f:
                c = csv.reader(f, delimiter='\t')
                #not_recognised = []
                for row in c:
                    #cid = row[0].split(':')
                    if not row[0][0]=='#':
                        #collect the info
                        mnxc = row[1]
                        name = mnxc_name[mnxc]
                        if len(row[0].split(':'))==1:
                            dbName = 'mnx'
                            dbCompId = row[0]
                        else:
                            dbName = row[0].split(':')[0]
                            dbCompId = ''.join(row[0].split(':')[1:])
                        #create the dicts
                        if not name in name_pubDB_xref:
                            name_pubDB_xref[name] = {}
                        if not dbName in name_pubDB_xref[name]:
                            name_pubDB_xref[name][dbName] = []
                        if not dbCompId in name_pubDB_xref[name][dbName]:
                            name_pubDB_xref[name][dbName].append(dbCompId)
        except FileNotFoundError:
            logging.error('compXref file not found')
            return {}
        return name_pubDB_xref


    ## Function to parse the chemp_prop.tsv file from MetanetX
    #
    #  Generate a dictionnary gaving the formula, smiles, inchi and inchikey for the components
    #
    #  @param self Object pointer
    #  @param chem_prop_path Input file path 
    #  @return mnxm_strc Dictionnary of formula, smiles, inchi and inchikey
    def mnxm_strc(self, chem_prop_path):
        #TODO: need to reduce the size of this file. As it stands its 250MB
        mnxm_strc = {}
        with open(chem_prop_path) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnxm_strc[row[0]] = {'forumla':  row[2], 'smiles': row[6], 'inchi': row[5], 'inchikey': row[8]}
                    #mnxm_strc[row[6]] = {'forumla':  row[2], 'mnxm': row[0], 'inchi': row[5], 'inchikey': row[8]}
                    for i in mnxm_strc[row[0]]:
                        if mnxm_strc[row[0]][i]=='' or mnxm_strc[row[0]][i]=='NA':
                            mnxm_strc[row[0]][i] = None
                    # if you have smiles
                    try:
                        if mnxm_strc[row[0]]['smiles']:
                            if not mnxm_strc[row[0]]['inchi']:
                                resConv = self._convert_depiction(idepic=mnxm_strc[row[0]]['smiles'], itype='smiles', otype={'inchi, inchikey'})
                                mnxm_strc[row[0]]['inchi'] = resConv['inchi']
                                mnxm_strc[row[0]]['inchikey'] = resConv['inchikey']
                        elif mnxm_strc[row[0]]['inchi']:
                            resConv = self._convert_depiction(idepic=mnxm_strc[row[0]]['inchi'], itype='inchi', otype={'smiles, inchikey'})
                            mnxm_strc[row[0]]['smiles'] = resConv['smiles']
                            mnxm_strc[row[0]]['inchikey'] = resConv['inchikey']
                    except (NotImplementedError, Exception) as e:
                        logging.warning('Could not convert the structures for '+str(row[0])+': '+str(mnxm_strc[row[0]]))
        return mnxm_strc


    ## Function exctract the dG of components
    #
    #
    #
    #  @param self Object pointer
    #  @param self.deprecatedMNXM_mnxm Dictionnary of old/new MNX identifiers
    #  @param chem_xref_path chem_xref.tsv file path
    #  @param cc_compounds_path cc_compounds.json.gz file path
    #  @param alberty_path alberty.json file path
    #  @param compounds_path compounds.csv file path
    def kegg_dG(self,
                cc_compounds_path,
                alberty_path,
                compounds_path):
        cc_alberty = {}
        ########################## compounds ##################
        #contains the p_kas and molecule decomposition 
        cid_comp = {}
        with open(compounds_path) as f:
            c = csv.reader(f, delimiter=',', quotechar='"')
            next(c)
            for row in c:
                cid_comp[row[-1].split(':')[1]] = {}
                cid_comp[row[-1].split(':')[1]]['atom_bag'] = literal_eval(row[3])
                cid_comp[row[-1].split(':')[1]]['p_kas'] = literal_eval(row[4])
                cid_comp[row[-1].split(':')[1]]['major_ms'] = int(literal_eval(row[6]))
                cid_comp[row[-1].split(':')[1]]['number_of_protons'] = literal_eval(row[7]) 
                cid_comp[row[-1].split(':')[1]]['charges'] = literal_eval(row[8])
        '''
        ###################### mnxm_kegg ################
        kegg_mnxm = {}
        with open(self._checkFilePath(chem_xref_path, 'chem_xref.tsv')) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                mnx = row[0].split(':') 
                if mnx[0]=='kegg' and mnx[1][0]=='C':
                    if mnx[1] in kegg_mnxm:
                        logging.warning(
                            'Replaced '+str(mnx[1])+': '+str(kegg_mnxm[mnx[1]])+' from '+str(row[1])
                        )
                    else:
                        try:
                            mnxm = self.deprecatedMNXM_mnxm[mnx[1]]
                        except KeyError:
                            mnxm = mnx[1]
                        kegg_mnxm[mnxm] = row[1]
        '''
        ####################### cc_compounds ############
        #TODO: seems like the new version of equilibrator got rid of this file... need to update the function
        #to take as input the new file --> i.e. the JSON input
        #notFound_cc = []
        gz_file = gzip.open(cc_compounds_path, 'rb')
        f_c = gz_file.read()
        c = json.loads(f_c)
        for cd in c:
            '''
            #find MNXM from CID
            try:
                mnxm = kegg_mnxm[cd['CID']]
            except KeyError:
                try:
                    mnxm = self.curated_kegg_mnxm[cd['CID']]
                    try:
                        mnxm = self.deprecatedMNXM_mnxm[mnxm]
                    except KeyError:
                        pass
                except KeyError:
                    logging.warning('Cannot find: '+str(cd))
                    notFound_cc.append(cd['CID'])
                    continue
            '''
            #find the compound descriptions
            try:
                cd.update(cid_comp[cd['CID']])
            except KeyError:
                pass
            #add the CID
            #if not mnxm in cc_alberty:
            if not cd['CID'] in cc_alberty:
                cc_alberty[cd['CID']] = {}
            if not 'component_contribution' in cc_alberty[cd['CID']]:
                cc_alberty[cd['CID']]['component_contribution'] = [cd]
            else:
                cc_alberty[cd['CID']]['component_contribution'].append(cd)
        ######################## alberty ################
        with open(alberty_path) as json_data:
            d = json.loads(json_data.read())
            for cd in d:
                '''
                #find the MNXM from CID
                try:
                    mnxm = kegg_mnxm[cd['cid']]
                except KeyError:
                    try:
                        mnxm = self.curated_kegg_mnxm[cd['cid']]
                        try:
                            mnxm = self.deprecatedMNXM_mnxm[mnxm]
                        except KeyError:
                            pass
                    except KeyError:
                        logging.warning('Cannot find: '+str(cd))
                        notFound_alberty.append(cd['cid'])
                        continue
                '''
                #find the compound description
                try:
                    cd.update(cid_comp[cd['CID']])
                except KeyError:
                    pass
                #add the CID
                #if not mnxm in cc_alberty:
                if not cd['CID'] in cc_alberty:
                    cc_alberty[cd['CID']] = {}
                if not 'alberty' in cc_alberty[cd['CID']]:
                    cc_alberty[cd['CID']]['alberty'] = [cd]
                else:
                    cc_alberty[cd['CID']]['alberty'].append(cd)
        return cc_alberty


    ## Function to parse the rules_rall.tsv file
    #
    #  Extract from the reactions rules the ruleID, the reactionID, the direction of the rule directed to the origin reaction
    #
    #  @param self The object pointer.
    #  @param path The input file path.
    #  @return rule Dictionnary describing each reaction rule
    def retro_reactions(self, path):
        try:
            with open(path, 'r') as f:
                reader = csv.reader(f, delimiter = '\t')
                next(reader)
                rule = {}
                for row in reader:
                    #rule[row[0]]= {'rule_id': row[0], 'rule_score': row[11], 'reaction':row[1], 'rel_direction': row[13], 'left': row[5], 'right': row[7]}
                    #NOTE: as of now all the rules are generated using MNX
                    #but it may be that other db are used, we are handling this case
                    #WARNING: can have multiple products so need to seperate them
                    products = {}
                    for i in row[7].split('.'):
                        if not i in products:
                            products[i] = 1
                        else:
                            products[i] += 1
                    try:
                        rule[row[0]] = {'rule_id': row[0], 'rule_score': float(row[11]), 'reac_id': row[1], 'subs_id': row[5], 'rel_direction': int(row[13]), 'left': {row[5]: 1}, 'right': products}
                    except ValueError:
                        logging.error('Problem converting rel_direction: '+str(row[13]))
                        logging.error('Problem converting rule_score: '+str(row[11]))
        except FileNotFoundError as e:
                logging.error('Could not read the rules_rall file ('+str(path)+')')
                return {}
        return rule


    def full_reac(self, path):
        #### for character matching that are returned
        DEFAULT_STOICHIO_RESCUE = {"4n": 4, "3n": 3, "2n": 2, 'n': 1,
                           '(n)': 1, '(N)': 1, '(2n)': 2, '(x)': 1,
                           'N': 1, 'm': 1, 'q': 1,
                           '0.01': 1, '0.1': 1, '0.5': 1, '1.5': 1,
                           '0.02': 1, '0.2': 1,
                           '(n-1)': 0, '(n-2)': -1}
        reaction = {}
        try:
            with open(path) as f:
                reader = csv.reader(f, delimiter='\t')
                next(reader)
                for row in reader:
                    tmp = {} # makes sure that if theres an error its not added
                    #parse the reaction equation
                    if len(row[1].split('='))==2:
                        #reac_left = row[1].split('=')[0]
                        #reac_right = row[1].split('=')[1]
                        ######### LEFT ######
                        #### MNX id
                        tmp['left'] = {}
                        for spe in re.findall('(\(n-1\)|\d+|4n|3n|2n|n|\(n\)|\(N\)|\(2n\)|\(x\)|N|m|q|\(n\-2\)|\d+\.\d+) ([\w\d]+)@\w+', row[1].split('=')[0]):
                            #1) try to rescue if its one of the values
                            try:
                                tmp['left'][spe[1]] = DEFAULT_STOICHIO_RESCUE[spe[0]]
                            except KeyError:
                                #2) try to convert to int if its not
                                try:
                                    tmp['left'][spe[1]] = int(spe[0])
                                except ValueError:
                                    logging.warning('Cannot convert '+str(spe[0]))
                                    continue
                        '''# Common name of chemicals -- Perhaps implement some day
                        #### chem names
                        tmp['chem_left'] = {}
                        for spe in re.findall('(\(n-1\)|\d+|4n|3n|2n|n|\(n\)|\(N\)|\(2n\)|\(x\)|N|m|q|\(n\-2\)|\d+\.\d+) `([^`]+)`', row[2].split('=')[0]):
                            #1) try to rescue if its one of the values
                            try:
                                tmp['chem_left'][spe[1]] = DEFAULT_STOICHIO_RESCUE[spe[0]]
                            except KeyError:
                                #2) try to convert to int if its not
                                try: 
                                    tmp['chem_left'][spe[1]] = int(spe[0])
                                except ValueError:
                                    logging.warning('Cannot convert '+str(spe[0]))
                                    continue
                        '''
                        ####### RIGHT #####
                        ####  MNX id
                        tmp['right'] = {}
                        for spe in re.findall('(\(n-1\)|\d+|4n|3n|2n|n|\(n\)|\(N\)|\(2n\)|\(x\)|N|m|q|\(n\-2\)|\d+\.\d+) ([\w\d]+)@\w+', row[1].split('=')[1]):
                            #1) try to rescue if its one of the values
                            try:
                                tmp['right'][spe[1]] = DEFAULT_STOICHIO_RESCUE[spe[0]]
                            except KeyError:
                                #2) try to convert to int if its not
                                try: 
                                    tmp['right'][spe[1]] = int(spe[0])
                                except ValueError:
                                    logging.warning('Cannot convert '+str(spe[0]))
                                    continue
                        ''' Common name of chemicals -- Perhaps implement some day
                        #### chem names
                        tmp['chem_right'] = {}
                        for spe in re.findall('(\(n-1\)|\d+|4n|3n|2n|n|\(n\)|\(N\)|\(2n\)|\(x\)|N|m|q|\(n\-2\)|\d+\.\d+) `([^`]+)`', row[2].split('=')[1]):
                            #1) try to rescue if its one of the values
                            try:
                                tmp['chem_right'][spe[1]] = DEFAULT_STOICHIO_RESCUE[spe[0]]
                            except KeyError:
                                #2) try to convert to int if its not
                                try: 
                                    tmp['chem_right'][spe[1]] = int(spe[0])
                                except ValueError:
                                    logging.warning('Cannot convert '+str(spe[0]))
                                    continue
                        '''
                        ####### DIRECTION ######
                        try:
                            tmp['direction'] = int(row[3])
                        except ValueError:
                            logging.error('Cannot convert '+str(row[3])+' to int')
                            continue
                        ### add the others
                        tmp['main_left'] = row[9].split(',')
                        tmp['main_right'] = row[10].split(',')
                        reaction[row[0]] = tmp
                    else:
                        logging.warning('There should never be more or less than a left and right of an euation')
                        continue
            return reaction
        except FileNotFoundError:
            logging.error('Cannot find file: '+str(path))
            return False
 


## Run all the functions
#
#  Run all the files required to generate the cache. Requires : mvc.db, chem_xref.tsv, chem_prop.tsv, cc_compounds.json.gz and alberty.json
#
#  @param self Object pointer
#TODO: change the input_cache to a compressed file with a given checksum to assure that it is fine
#TODO: consider checksumming the individual files that are generated from here
if __name__ == "__main__":
    if not os.path.isdir(os.getcwd()+'/cache'):
        os.mkdir('cache')
    #dirname = os.path.dirname(os.path.abspath( __file__ ))
    cache = rpCache()
    #cache.deprecatedMNXM_mnxm
    logging.info('Generating deprecatedMNXM_mnxm')
    cache.deprecatedMNXM_mnxm = cache.deprecatedMNXM('input_cache/chem_xref.tsv')
    pickle.dump(cache.deprecatedMNXM_mnxm, open('cache/deprecatedMNXM_mnxm.pickle', 'wb'))
    #mnxm_dG
    logging.info('Generating mnxm_dG')
    pickle.dump(cache.kegg_dG('input_cache/cc_compounds.json.gz',
        'input_cache/alberty.json',
        'input_cache/compounds.csv'),
        open('cache/kegg_dG.pickle', 'wb'))
    #rr_reactions
    logging.info('Generating rr_reactions')
    rr_reactions = cache.retro_reactions('input_cache/rules_rall.tsv')
    pickle.dump(rr_reactions, open('cache/rr_reactions.pickle', 'wb'))
    #full_reactions
    logging.info('Generating full_reactions')
    pickle.dump(cache.full_reac('input_cache/rxn_recipes.tsv'), 
            open('cache/full_reactions.pickle', 'wb'))
    #mnxm_strc --> use gzip since it is a large file
    logging.info('Parsing the SMILES and InChI')
    #pickle.dump(cache.mnxm_strc(), open('cache/mnxm_strc.pickle', 'wb'))
    pickle.dump(cache.mnxm_strc('input_cache/chem_prop.tsv'), 
            gzip.open('cache/mnxm_strc.pickle.gz','wb'))
    #xref --> use gzip since it is a large file
    logging.info('Parsing the Cross-references')
    pickle.dump(cache.chemXref('input_cache/chem_xref.tsv'), gzip.open('cache/chemXref.pickle.gz','wb'))
    pickle.dump(cache.reacXref('input_cache/reac_xref.tsv', 'input_cache/reac_prop.tsv'), gzip.open('cache/reacXref.pickle.gz','wb'))
    pickle.dump(cache.compXref('input_cache/comp_xref.tsv', 'input_cache/comp_prop.tsv'), gzip.open('cache/compXref.pickle.gz','wb'))
    #copy the other file required
    copyfile('input_cache/cc_preprocess.npz', 'cache/cc_preprocess.npz')

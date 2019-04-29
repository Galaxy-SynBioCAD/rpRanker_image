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
    def __init__(self, inputPath=None):
        if inputPath and inputPath[-1:]=='/':
            inputPath = inputPath[:-1]
        self.globalPath = inputPath
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

    #######################################################################
    ############################# PRIVATE FUNCTIONS ####################### 
    #######################################################################

    def _checkFilePath(self, path, filename):
        """Check that the directory and the filename are valid and choose to use
        either the local or the global path
        """
        if path==None:
            if self.globalPath==None:
                logging.error('Both global path and local are not set')
                return None
            else:
                if os.path.isdir(self.globalPath):
                    try:
                        fName = [i for i in os.listdir(self.globalPath) 
                            if not i.find(filename)==-1 
                            if not i[-3:]=='swp'
                            if not i[-1]=='#'][0]
                        logging.info('Automatically selected '+str(fName))
                    except IndexError:
                        logging.error('Problem finding the correct '+str(filename)+' in '+str(path))
                        return None
                    if os.path.isfile(self.globalPath+'/'+fName):
                        return self.globalPath+'/'+fName
                    else:
                        logging.error('Global path file: '+str(fName)+', does not exist')
                        return None
                else:
                    logging.error('Global path is not a directory: '+str(self.globalPath))
                    return None
        else:
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

    ########################################################
    ############################# FUNCTIONS ################
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
    def deprecatedMNXM(self, chem_xref_path=None):
        a = {}
        with open(self._checkFilePath(chem_xref_path, 'chem_xref.tsv')) as f:
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
    def chem_xref(self, chem_xref_path=None):
        chemXref = {}
        with open(self._checkFilePath(chem_xref_path, 'chem_xref.tsv')) as f:
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


    ## Function to parse the reac_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of old to new MetanetX identifiers
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return a The dictionnary of identifiers  
    #TODO: save the self.deprecatedMNXM_mnxm to be used in case there rp_paths uses an old version of MNX
    def reac_xref(self, reac_xref_path=None):
        reacXref = {}
        with open(self._checkFilePath(reac_xref_path, 'reac_xref.tsv')) as f:
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
        return reacXref


    ## Function to parse the comp_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of old to new MetanetX identifiers
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return a The dictionnary of identifiers  
    #TODO: save the self.deprecatedMNXM_mnxm to be used in case there rp_paths uses an old version of MNX
    def comp_xref(self, comp_xref_path=None, comp_prop_path=None):
        mnxc_name = {}
        with open(self._checkFilePath(comp_prop_path, 'comp_prop.tsv')) as f:
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
            with open(self._checkFilePath(comp_xref_path, 'comp_xref.tsv')) as f:
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
            logging.error('comp_xref file not found')
            return {}
        return name_pubDB_xref

    ## Function to parse the chemp_prop.tsv file from MetanetX
    #
    #  Generate a dictionnary gaving the formula, smiles, inchi and inchikey for the components
    #
    #  @param self Object pointer
    #  @param chem_prop_path Input file path 
    #  @return smiles_inchi Dictionnary of formula, smiles, inchi and inchikey
    def smiles_inchi(self, chem_prop_path=None):
        #TODO: need to reduce the size of this file. As it stands its 250MB
        smiles_inchi = {}
        with open(self._checkFilePath(chem_prop_path, 'chem_prop.tsv')) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    smiles_inchi[row[0]] = {'forumla':  row[2], 'smiles': row[6], 'inchi': row[5], 'inchikey': row[8]}
        return smiles_inchi

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
                cc_compounds_path=None,
                alberty_path=None,
                compounds_path=None):
        cc_alberty = {}
        ########################## compounds ##################
        #contains the p_kas and molecule decomposition 
        cid_comp = {}
        with open(self._checkFilePath(compounds_path, 'compounds.csv')) as f:
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
        gz_file = gzip.open(self._checkFilePath(cc_compounds_path, 'cc_compounds.json.gz'), 'rb')
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
        with open(self._checkFilePath(alberty_path, 'alberty.json')) as json_data:
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
    def retro_reactions(self, path=None):
        try:
            with open(self._checkFilePath(path, 'rules_rall'), 'r') as f:
                reader = csv.reader(f, delimiter = '\t')
                next(reader)
                rule = {}
                for row in reader:
                    rule[row[0]]= {'reaction':row[1], 'rel_direction': row[13], 'left': row[5], 'right': row[7]}
        except (TypeError, FileNotFoundError) as e:
                logging.error('Could not read the rules_rall file ('+str(path)+')')
                return {}
        return(rule)


    ## Function to parse rxn_recipes.tsv file
    #
    #  Extract the substracts and products of the origin reaction with the stochiometry
    #
    #  @param self The object pointer.
    #  @param self.deprecatedMNXM_mnxm Dictionnary of old to new version of MetanetX identifiers
    #  @param path The input file path.
    #  @return reaction Dictionnnary containing the description of the reaction of origin for each reactionID
    def full_reac(self, path=None):
        try:
            with open(self._checkFilePath(path, 'rxn_recipes')) as f:
                def dico(liste, i=0, x=1):
                        while i <= len(liste)-1:
                            tmp = liste[i];
                            liste[i] = liste[x]
                            liste[x] = tmp
                            i = i+2
                            x = x+2
                        dico = dict(itertools.zip_longest(*[iter(liste)] * 2, fillvalue=""))
                        return dico
                reader = csv.reader(f, delimiter = '\t')
                next(reader)
                reaction = {}
                for row in reader:
                    reaction[row[0]]= {'main_left': row[9], 'main_right': row[10], 'left':{}, 'right':{}, 'direction': row[3]}
                    #tmp_l = (((((((((row[1].split('=')[0]).replace('+ ',"")).replace('@MNXD1','')).replace('@MNXD2','')).replace('@MNXDX','')).replace('@BOUNDARY','')).replace('@\S',''))).split(' '))[:-1]
                    tmp_l = ((row[1].split('=')[0]).replace('+ ',"")).split(' ')[:-1]
                    #tmp_r = ((((((((row[1].split('=')[1]).replace('+ ',"")).replace('@MNXD1','')).replace('@MNXD2','')).replace('@MNXDX','')).replace('@BOUNDARY',''))).split(' '))[1:]
                    tmp_r = ((row[1].split('=')[1]).replace('+ ',"")).split(' ')[1:]
                    
                    p = re.compile('@\\w+')
                    for x,i in enumerate(tmp_l):
                        l = re.sub(p, '', i)
                        tmp_l[x]= l
                    for x,i in enumerate(tmp_r):
                        r = re.sub(p, '', i)
                        tmp_r[x]= r
                    
                    dico_l = dico(tmp_l)
                    dico_r = dico(tmp_r)
                    
                    #remove the main_right and the main_left because already in the out_path.csv
                    #WARNING: There is a small chance that you will remove elements with MNX codes that 
                    #are in fact not described as the CMP
                    '''for i in reaction[row[0]]['main_right'].split(','):
                        del dico_r[i]
                    for i in reaction[row[0]]['main_left'].split(','):
                        del dico_l[i]'''

                    for i in dico_l:
                        tmp_i = i
                        while tmp_i in self.deprecatedMNXM_mnxm:
                            tmp_i = self.deprecatedMNXM_mnxm[tmp_i]
                        reaction[row[0]]['left'][tmp_i] = dico_l[i]
                    for i in dico_r:
                        tmp_i = i
                        while tmp_i in self.deprecatedMNXM_mnxm:
                            tmp_i = self.deprecatedMNXM_mnxm[tmp_i]
                        reaction[row[0]]['right'][tmp_i] = dico_r[i]
        except (TypeError, FileNotFoundError) as e:
                logging.error('Could not read the rxn_recipes file ('+str(path)+')')
                return {}
        return reaction



## Run all the functions
#
#  Run all the files required to generate the cache. Requires : mvc.db, chem_xref.tsv, chem_prop.tsv, cc_compounds.json.gz and alberty.json
#
#  @param self Object pointer
#TODO: change the input_cache to a compressed file with a given checksum to assure that it is fine
#TODO: consider checksumming the individual files that are generated from here
if __name__ == "__main__":
    cache = rpCache(os.path.abspath('input_cache'))
    if cache.globalPath==None:
        logging.error('Need to define the global path to use all()')
        sys.exit()
    if not os.path.isdir(os.getcwd()+'/cache'):
        os.mkdir('cache')
    #cache.deprecatedMNXM_mnxm
    logging.info('Generating deprecatedMNXM_mnxm')
    cache.deprecatedMNXM_mnxm = cache.deprecatedMNXM()
    pickle.dump(cache.deprecatedMNXM_mnxm, open('cache/deprecatedMNXM_mnxm.pickle', 'wb'))
    #mnxm_dG
    logging.info('Generating mnxm_dG')
    pickle.dump(cache.kegg_dG(), open('cache/kegg_dG.pickle', 'wb'))
    #rr_reactions
    logging.info('Generating rr_reactions')
    rr_reactions = cache.retro_reactions()
    pickle.dump(rr_reactions, open('cache/rr_reactions.pickle', 'wb'))
    #full_reactions
    logging.info('Generating full_reactions')
    pickle.dump(cache.full_reac(), open('cache/full_reactions.pickle', 'wb'))
    #smiles_inchi --> use gzip since it is a large file
    logging.info('Parsing the SMILES and InChI')
    #pickle.dump(cache.smiles_inchi(), open('cache/smiles_inchi.pickle', 'wb'))
    pickle.dump(cache.smiles_inchi(), gzip.open('cache/smiles_inchi.pickle.gz','wb'))
    #xref --> use gzip since it is a large file
    logging.info('Parsing the Cross-references')
    pickle.dump(cache.chem_xref(), gzip.open('cache/chemXref.pickle.gz','wb'))
    pickle.dump(cache.reac_xref(), gzip.open('cache/reacXref.pickle.gz','wb'))
    pickle.dump(cache.comp_xref(), gzip.open('cache/compXref.pickle.gz','wb'))
    '''
    pub_mnx_chem_xref, mnx_pub_chem_xref = cache.chem_xref()
    pickle.dump(pub_mnx_chem_xref, gzip.open('cache/pub_mnx_chem_xref.gz','wb'))
    pickle.dump(mnx_pub_chem_xref, gzip.open('cache/mnx_pub_chem_xref.gz','wb'))
    pub_mnx_reax_xref, mnx_pub_reac_xref = cache.reac_xref()
    pickle.dump(pub_mnx_reax_xref, gzip.open('cache/pub_mnx_reax_xref.gz','wb'))
    pickle.dump(mnx_pub_reac_xref, gzip.open('cache/mnx_pub_reac_xref.gz','wb'))
    pub_mnx_comp_xref, mnx_pub_comp_xref = cache.comp_xref()
    pickle.dump(pub_mnx_comp_xref, gzip.open('cache/pub_mnx_comp_xref.gz','wb'))
    pickle.dump(mnx_pub_comp_xref, gzip.open('cache/mnx_pub_comp_xref.gz','wb'))
    '''
    #pickle.dump(cache.xref(), gzip.open('cache/Id_xref.pickle.gz','wb'))

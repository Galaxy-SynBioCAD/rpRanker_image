import csv
import gzip
import os
import json
import pickle
#import sqlite3
import re
import itertools
from ast import literal_eval
import gzip
from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi, MolToInchiKey, AddHs
from shutil import copyfile
import tarfile

import logging

## @package rpCache
#
# These are a collection of functions that parse description files of reactions (MNX), chemical species (MNX) and model compartments (MNX) to be used by rpRanker. Furthermore, it parses euilibrator files to calculate the thermodynamics of reactions. This includes pre-calculated values as well as on the fly. Lastly reaction rule files are parsed to reconstruct full reactions from monocomponent reactions.


#######################################################
################### USER DEFINED ERROR ################
#######################################################


## Error function for the convertion of structures
#
class Error(Exception):
    pass


## Error function for the convertion of structures
#
class DepictionError(Error):
    def __init__(self, message):
        #self.expression = expression
        self.message = message


## Class to generate the cache
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
        self.logger = logging.getLogger(__name__)
        self.logger.info('Started instance of rpCache')
        self.convertMNXM = {'MNXM162231': 'MNXM6',
                'MNXM84': 'MNXM15',
                'MNXM96410': 'MNXM14',
                'MNXM114062': 'MNXM3',
                'MNXM145523': 'MNXM57',
                'MNXM57425': 'MNXM9',
                'MNXM137': 'MNXM588022'}
        self.deprecatedMNXM_mnxm = {}
        self.deprecatedMNXR_mnxr = {}


    #######################################################
    ################### PRIVATE FUNCTION ##################
    #######################################################


    ## Convert chemical depiction to others type of depictions
    #
    # Usage example:
    # - convert_depiction(idepic='CCO', otype={'inchi', 'smiles', 'inchikey'})
    # - convert_depiction(idepic='InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', itype='inchi', otype={'inchi', 'smiles', 'inchikey'})
    #  @param self The object pointer
    #  @param idepic String depiction to be converted, str
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
            raise DepictionError('Import error from depiction "{}" of type "{}"'.format(idepic, itype))
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


    #[TODO] merge the two functions
    ## Function to parse the chem_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of old to new MetanetX identifiers to make sure that we always use the freshest id's.
    # This can include more than one old id per new one and thus returns a dictionnary. Private function
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return Dictionnary of identifiers  
    #TODO: save the self.deprecatedMNXM_mnxm to be used in case there rp_paths uses an old version of MNX
    def _deprecatedMNXM(self, chem_xref_path):
        self.deprecatedMNXM_mnxm = {}
        with open(chem_xref_path) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnx = row[0].split(':')
                    if mnx[0]=='deprecated':
                        self.deprecatedMNXM_mnxm[mnx[1]] = row[1]
            self.deprecatedMNXM_mnxm.update(self.convertMNXM)
            self.deprecatedMNXM_mnxm['MNXM01'] = 'MNXM1'


    ## Function to parse the reac_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of old to new MetanetX identifiers to make sure that we always use the freshest id's.
    # This can include more than one old id per new one and thus returns a dictionnary. Private function
    #
    #  @param self Object pointer
    #  @param reac_xref_path Input file path
    #  @return Dictionnary of identifiers  
    def _deprecatedMNXR(self, reac_xref_path):
        self.deprecatedMNXMR_mnxr = {}
        with open(reac_xref_path) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnx = row[0].split(':')
                    if mnx[0]=='deprecated':
                        self.deprecatedMNXR_mnxr[mnx[1]] = row[1]


    ## Function to create a dictionnary of old to new chemical id's 
    #
    #  Generate a one-to-one dictionnary of old id's to new ones. Private function
    #
    # TODO: check other things about the mnxm emtry like if it has the right structure etc...
    def _checkMNXMdeprecated(self, mnxm):
        try:
            return self.deprecatedMNXM_mnxm[mnxm]
        except KeyError:
            return mnxm


    ## Function to create a dictionnary of old to new reaction id's
    #
    # TODO: check other things about the mnxm emtry like if it has the right structure etc...
    def _checkMNXRdeprecated(self, mnxr):
        try:
            return self.deprecatedMNXR_mnxr[mnxr]
        except KeyError:
            return mnxr


    #NOTE: no need to do the same for compartments :)


    ########################################################
    ####################### PUBLIC FUNCTIONS ###############
    ######################################################## 


    ################### MetaNetX ############################   


    ## Function to parse the chemp_prop.tsv file from MetanetX and compounds.tsv from RetroRules. Uses the InchIkey as key to the dictionnary
    #
    #  Generate a dictionnary gaving the formula, smiles, inchi and inchikey for the components
    #
    #  @param self Object pointer
    #  @param chem_prop_path Input file path 
    #  @return mnxm_strc Dictionnary of formula, smiles, inchi and inchikey
    def mnx_strc(self, rr_compounds_path, chem_prop_path):
        mnxm_strc = {}
        inchikey_mnxm = {}
        #with open(rr_compounds_path) as f:
        #    c = csv.reader(f, delimiter='\t')
        #    for row in c:
        for row in csv.DictReader(open(rr_compounds_path), delimiter='\t'):
            #if not row[0][0]=='#':
            tmp = {'forumla':  None, 
                    'smiles': None, 
                    'inchi': row['inchi'], 
                    'inchikey': None, 
                    'mnxm': self._checkMNXMdeprecated(row['cid']), 
                    'name': None}
            try:
                resConv = self._convert_depiction(idepic=tmp['inchi'], itype='inchi', otype={'smiles','inchikey'})
                for i in resConv:
                    tmp[i] = resConv[i]
            except DepictionError as e:
                self.logger.warning('Could not convert some of the structures: '+str(tmp))
                self.logger.warning(e)
            mnxm_strc[tmp['mnxm']] = tmp
        with open(chem_prop_path) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnxm = self._checkMNXMdeprecated(row[0])
                    tmp = {'forumla':  row[2], 
                            'smiles': row[6], 
                            'inchi': row[5], 
                            'inchikey': row[8], 
                            'mnxm': mnxm, 
                            'name': row[1]}
                    for i in tmp:
                        if tmp[i]=='' or tmp[i]=='NA':
                            tmp[i] = None
                    if mnxm in mnxm_strc:
                        mnxm_strc[mnxm]['forumla'] = row[2]
                        mnxm_strc[mnxm]['name'] = row[1]
                        if not mnxm_strc[mnxm]['smiles'] and tmp['smiles']:
                            mnxm_strc[mnxm]['smiles'] = tmp['smiles']
                        if not mnxm_strc[mnxm]['inchikey'] and tmp['inchikey']:
                            mnxm_strc[mnxm]['inchikey'] = tmp['inchikey']
                        #### inchikey_mnxm ###
                        if not tmp['inchikey'] in inchikey_mnxm:
                            inchikey_mnxm[tmp['inchikey']] = {'mnx': []}
                        inchikey_mnxm[tmp['inchikey']]['mnx'].append(tmp['mnxm'])
                    else:
                        #check to see if the inchikey is valid or not
                        otype = set({})
                        if not tmp['inchikey']:
                            otype.add('inchikey')
                        if not tmp['smiles']:
                            otype.add('smiles')
                        if not tmp['inchi']:
                            otype.add('inchi')
                        itype = ''
                        if tmp['inchi']:
                            itype = 'inchi'
                        elif tmp['smiles']:
                            itype = 'smiles'
                        else:
                            self.logger.warning('No valid entry for the convert_depiction function')
                            continue
                        try:
                            resConv = self._convert_depiction(idepic=tmp[itype], itype=itype, otype=otype)
                            for i in resConv:
                                tmp[i] = resConv[i]
                        except DepictionError as e:
                            self.logger.warning('Could not convert some of the structures: '+str(tmp))
                            self.logger.warning(e)
                        mnxm_strc[tmp['mnxm']] = tmp
                        #### inchikey_mnxm ###
                        if not tmp['inchikey'] in inchikey_mnxm:
                            inchikey_mnxm[tmp['inchikey']] = {'mnx': []}
                        inchikey_mnxm[tmp['inchikey']]['mnx'].append(tmp['mnxm'])
        return mnxm_strc, inchikey_mnxm


    ## Function to parse the chem_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of all cross references for a given chemical id (MNX) to other database id's
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return a The dictionnary of identifiers  
    #TODO: save the self.deprecatedMNXM_mnxm to be used in case there rp_paths uses an old version of MNX
    def mnx_chemXref(self, chem_xref_path):
        chemXref = {}
        with open(chem_xref_path) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnx = self._checkMNXMdeprecated(row[1])
                    if len(row[0].split(':'))==1:
                        dbName = 'mnx'
                        dbId = row[0]
                    else:
                        dbName = row[0].split(':')[0]
                        dbId = ''.join(row[0].split(':')[1:])
                        if dbName=='deprecated':
                            dbName = 'mnx'
                    #mnx
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


    ## Function to parse the reacXref.tsv file of MetanetX and rxn_recipes.tsv from RetroRules
    #
    #  Generate a dictionnary of all cross references for a given reaction id (MNX) to other database id's 
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @param rxn_recipes_path Input file path for the reaction recipes
    #  @return a The dictionnary of identifiers  
    def mnx_reacXref(self, reacXref_path, rxn_recipes_path):
        reacXref = {}
        with open(reacXref_path) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnx = self._checkMNXRdeprecated(row[1])
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
        #with open(rxn_recipes_path) as f:
        #    c = csv.reader(f, delimiter='\t')
        #    for row in c:
        for row in csv.DictReader(open(rxn_recipes_path), delimiter='\t'):
            if row['EC_number']=='':
                #### WARNING: this is temporary TODO TODO TODO TODO once updated rxn_recipes
                #row[0] as MNX that is
                mnx = self._checkMNXRdeprecated(row['#Reaction_ID'])
                if not mnx in reacXref:
                    reacXref[mnx] = {}
                dbName = 'ec'
                if not dbName in reacXref[mnx]:
                    reacXref[mnx][dbName] = []
                for ec in row['EC_number'].split(';'):
                    if not ec in reacXref[mnx][dbName]:
                        reacXref[mnx][dbName].append(ec)
        return reacXref


    ## Function to parse the compXref.tsv file of MetanetX
    #
    #  Generate a dictionnary of compartments id's (MNX) to other database id's
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return a The dictionnary of identifiers  
    #TODO: save the self.deprecatedMNXM_mnxm to be used in case there rp_paths uses an old version of MNX
    def mnx_compXref(self, compXref_path):
        name_pubDB_xref = {}
        compName_mnxc = {}
        try:
            with open(compXref_path) as f:
                c = csv.reader(f, delimiter='\t')
                #not_recognised = []
                for row in c:
                    #cid = row[0].split(':')
                    if not row[0][0]=='#':
                        #collect the info
                        mnxc = row[1]
                        if len(row[0].split(':'))==1:
                            dbName = 'mnx'
                            dbCompId = row[0]
                        else:
                            dbName = row[0].split(':')[0]
                            dbCompId = ''.join(row[0].split(':')[1:])
                            dbCompId = dbCompId.lower()
                        if dbName=='deprecated':
                            dbName = 'mnx'
                        #create the dicts
                        if not mnxc in name_pubDB_xref:
                            name_pubDB_xref[mnxc] = {}
                        if not dbName in name_pubDB_xref[mnxc]:
                            name_pubDB_xref[mnxc][dbName] = []
                        if not dbCompId in name_pubDB_xref[mnxc][dbName]:
                            name_pubDB_xref[mnxc][dbName].append(dbCompId)
                        #create the reverse dict
                        if not dbCompId in compName_mnxc:
                            compName_mnxc[dbCompId] = mnxc
        except FileNotFoundError:
            self.logger.error('compXref file not found')
            return {}
        return name_pubDB_xref, compName_mnxc


    ####################### Equilibrator ###############################


    ## Function exctract the dG of components
    #
    #  This function parses a file from component analysis and uses KEGG id's. This means that to use precalculated
    # values in this file a given ID must have a KEGG id listed here
    #
    #  @param self Object pointer
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
                        self.logger.warning(
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
                    self.logger.warning('Cannot find: '+str(cd))
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
                        self.logger.warning('Cannot find: '+str(cd))
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

    
    ####################################### RetroRules #####################################


    ## Function to parse the rules_rall.tsv from RetroRules
    #
    #  Extract from the reactions rules the ruleID, the reactionID, the direction of the rule directed to the origin reaction
    #
    #  @param self The object pointer.
    #  @param path The input file path.
    #  @return rule Dictionnary describing each reaction rule
    def retro_reactions(self, rules_rall_path):
        try:
            #with open(rules_rall_path, 'r') as f:
            #    reader = csv.reader(f, delimiter = '\t')
            #    next(reader)
            #    rule = {}
            #    for row in reader:
            rule = {}
            for row in csv.DictReader(open(rules_rall_path), delimiter='\t'):
                #rule[row[0]]= {'rule_id': row[0], 'rule_score': row[11], 'reaction':row[1], 'rel_direction': row[13], 'left': row[5], 'right': row[7]}
                #NOTE: as of now all the rules are generated using MNX
                #but it may be that other db are used, we are handling this case
                #WARNING: can have multiple products so need to seperate them
                products = {}
                for i in row['Product_IDs'].split('.'):
                    mnxm = self._checkMNXMdeprecated(i)
                    if not mnxm in products:
                        products[mnxm] = 1
                    else:
                        products[mnxm] += 1
                try:
                    '''
                    rule[row['# Rule_ID']] = {'rule_id': row['# Rule_ID'], 'rule_score': float(row['Score_normalized']), 'reac_id': self._checkMNXRdeprecated(row['Reaction_ID']), 'subs_id': self._checkMNXMdeprecated(row['Substrate_ID']), 'rel_direction': int(row['Rule_relative_direction']), 'left': {self._checkMNXMdeprecated(row['Substrate_ID']): 1}, 'right': products}
                    '''
                    #WARNING: one reaction rule can have multiple reactions associated with them
                    #To change when you can set subpaths from the mutliple numbers of 
                    #we assume that the reaction rule has multiple unique reactions associated
                    if row['# Rule_ID'] not in rule:
                        rule[row['# Rule_ID']] = {}
                    if row['# Rule_ID'] in rule[row['# Rule_ID']]:
                        self.logger.warning('There is already reaction '+str(row['# Rule_ID'])+' in reaction rule '+str(row['# Rule_ID']))
                    rule[row['# Rule_ID']][row['Reaction_ID']] = {'rule_id': row['# Rule_ID'], 'rule_score': float(row['Score_normalized']), 'reac_id': self._checkMNXRdeprecated(row['Reaction_ID']), 'subs_id': self._checkMNXMdeprecated(row['Substrate_ID']), 'rel_direction': int(row['Rule_relative_direction']), 'left': {self._checkMNXMdeprecated(row['Substrate_ID']): 1}, 'right': products}
                except ValueError:
                    self.logger.error('Problem converting rel_direction: '+str(row['Rule_relative_direction']))
                    self.logger.error('Problem converting rule_score: '+str(row['Score_normalized']))
        except FileNotFoundError as e:
                self.logger.error('Could not read the rules_rall file ('+str(path)+')')
                return {}
        return rule

    
    ## Generate complete reactions from the rxn_recipes.tsv from RetroRules
    #
    #  These are the compplete reactions from which the reaction rules are generated from. This is used to
    # reconstruct the full reactions from monocomponent reactions
    #
    #  @param self The pointer object
    #  @param rxn_recipes_path Path to the recipes file
    #  @return Boolean that determines the success or failure of the function
    def full_reac(self, rxn_recipes_path):
        #### for character matching that are returned
        DEFAULT_STOICHIO_RESCUE = {"4n": 4, "3n": 3, "2n": 2, 'n': 1,
                           '(n)': 1, '(N)': 1, '(2n)': 2, '(x)': 1,
                           'N': 1, 'm': 1, 'q': 1,
                           '0.01': 1, '0.1': 1, '0.5': 1, '1.5': 1,
                           '0.02': 1, '0.2': 1,
                           '(n-1)': 0, '(n-2)': -1}
        reaction = {}
        try:
            #with open(rxn_recipes_path) as f:
            #    reader = csv.reader(f, delimiter='\t')
            #    next(reader)
            #    for row in reader:
            for row in csv.DictReader(open(rxn_recipes_path), delimiter='\t'):
                tmp = {} # makes sure that if theres an error its not added
                #parse the reaction equation
                if not len(row['Equation'].split('='))==2:
                    self.logger.warning('There should never be more or less than a left and right of an euation')
                    self.logger.warnin(row['Equation'])
                    continue 
                #reac_left = row[1].split('=')[0]
                #reac_right = row[1].split('=')[1]
                ######### LEFT ######
                #### MNX id
                tmp['left'] = {}
                for spe in re.findall('(\(n-1\)|\d+|4n|3n|2n|n|\(n\)|\(N\)|\(2n\)|\(x\)|N|m|q|\(n\-2\)|\d+\.\d+) ([\w\d]+)@\w+', row['Equation'].split('=')[0]):
                    #1) try to rescue if its one of the values
                    try:
                        tmp['left'][self._checkMNXMdeprecated(spe[1])] = DEFAULT_STOICHIO_RESCUE[spe[0]]
                    except KeyError:
                        #2) try to convert to int if its not
                        try:
                            tmp['left'][self._checkMNXMdeprecated(spe[1])] = int(spe[0])
                        except ValueError:
                            self.logger.warning('Cannot convert '+str(spe[0]))
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
                            self.logger.warning('Cannot convert '+str(spe[0]))
                            continue
                '''
                ####### RIGHT #####
                ####  MNX id
                tmp['right'] = {}
                for spe in re.findall('(\(n-1\)|\d+|4n|3n|2n|n|\(n\)|\(N\)|\(2n\)|\(x\)|N|m|q|\(n\-2\)|\d+\.\d+) ([\w\d]+)@\w+', row['Equation'].split('=')[1]):
                    #1) try to rescue if its one of the values
                    try:
                        tmp['right'][self._checkMNXMdeprecated(spe[1])] = DEFAULT_STOICHIO_RESCUE[spe[0]]
                    except KeyError:
                        #2) try to convert to int if its not
                        try: 
                            tmp['right'][self._checkMNXMdeprecated(spe[1])] = int(spe[0])
                        except ValueError:
                            self.logger.warning('Cannot convert '+str(spe[0]))
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
                            self.logger.warning('Cannot convert '+str(spe[0]))
                            continue
                '''
                ####### DIRECTION ######
                try:
                    tmp['direction'] = int(row['Direction'])
                except ValueError:
                    self.logger.error('Cannot convert '+str(row['Direction'])+' to int')
                    continue
                ### add the others
                tmp['main_left'] = row['Main_left'].split(',')
                tmp['main_right'] = row['Main_right'].split(',')
                reaction[self._checkMNXRdeprecated(row['#Reaction_ID'])] = tmp
            return reaction
        except FileNotFoundError:
            self.logger.error('Cannot find file: '+str(path))
            return False
 

## Run all the functions
#
#  Run all the files required to generate the cache. Requires : mvc.db, chem_xref.tsv, chem_prop.tsv, cc_compounds.json.gz and alberty.json
#
#TODO: change the input_cache to a compressed file with a given checksum to assure that it is fine
#TODO: consider checksumming the individual files that are generated from here
if __name__ == "__main__":
    if not os.path.isdir(os.getcwd()+'/cache'):
        os.mkdir('cache')
    #dirname = os.path.dirname(os.path.abspath( __file__ ))
    cache = rpCache()
    #generate the deprecated
    #self.logger.info('Generating deprecatedMNXM_mnxm')
    cache._deprecatedMNXM('input_cache/chem_xref.tsv')
    cache._deprecatedMNXR('input_cache/reac_xref.tsv')
    pickle.dump(cache.deprecatedMNXM_mnxm, open('cache/deprecatedMNXM_mnxm.pickle', 'wb'))
    pickle.dump(cache.deprecatedMNXR_mnxr, open('cache/deprecatedMNXR_mnxr.pickle', 'wb'))
    #mnxm_dG
    #self.logger.info('Generating mnxm_dG')
    pickle.dump(cache.kegg_dG('input_cache/cc_compounds.json.gz',
        'input_cache/alberty.json',
        'input_cache/compounds.csv'),
        open('cache/kegg_dG.pickle', 'wb'))
    #rr_reactions
    #self.logger.info('Generating rr_reactions')
    rr_reactions = cache.retro_reactions('input_cache/rules_rall.tsv')
    pickle.dump(rr_reactions, open('cache/rr_reactions.pickle', 'wb'))
    #full_reactions
    #self.logger.info('Generating full_reactions')
    pickle.dump(cache.full_reac('input_cache/rxn_recipes.tsv'), 
            open('cache/full_reactions.pickle', 'wb'))
    #mnxm_strc --> use gzip since it is a large file
    #self.logger.info('Parsing the SMILES and InChI')
    #pickle.dump(cache.mnxm_strc(), open('cache/mnxm_strc.pickle', 'wb'))
    mnx_strc, inchikey_mnxm = cache.mnx_strc('input_cache/compounds.tsv', 'input_cache/chem_prop.tsv')
    pickle.dump(mnx_strc, gzip.open('cache/mnxm_strc.pickle.gz','wb'))
    pickle.dump(inchikey_mnxm, gzip.open('cache/inchikey_mnxm.pickle.gz','wb'))
    #xref --> use gzip since it is a large file
    #self.logger.info('Parsing the Cross-references')
    pickle.dump(cache.mnx_chemXref('input_cache/chem_xref.tsv'), gzip.open('cache/chemXref.pickle.gz','wb'))
    pickle.dump(cache.mnx_reacXref('input_cache/reac_xref.tsv', 'input_cache/rxn_recipes.tsv'), gzip.open('cache/reacXref.pickle.gz','wb'))
    name_pubDB_xref, compName_mnxc = cache.mnx_compXref('input_cache/comp_xref.tsv')
    pickle.dump(name_pubDB_xref, gzip.open('cache/compXref.pickle.gz','wb'))
    pickle.dump(compName_mnxc, gzip.open('cache/nameCompXref.pickle.gz','wb'))
    #copy the other file required
    copyfile('input_cache/cc_preprocess.npz', 'cache/cc_preprocess.npz')
    #copy the rules_rall.tsv 
    #with tarfile.open('cache/rules_rall.tsv.tar.xz', "w:xz") as tar:
    #    tar.add('input_cache/rules_rall.tsv')


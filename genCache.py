import csv
import logging
import json
import gzip
import os
import pickle
import sqlite3

class Cache:
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
        #personally lookes like the KEGG to MNXM
        self.curated_kegg_mnxm = {'C80055': 'MNXM1029', 'C80033': 'MNXM4958', 'C80027': 'MNXM14271', 'C80026': 'MNXM37357', 'C80024': 'MNXM37357', 'C80021' :'MNXM34393', 'C18385': 'MNXM155504','C02747': 'MNXM3752', 'C00281': 'MNXM271', 'C00660': 'MNXM7242', 'C02229': 'MNXM493', 'C02799': 'MNXM148650', 'C03009': 'MNXM64052', 'C03293': 'MNXM4718', 'C03543': 'MNXM5043', 'C03857': 'MNXM48439', 'C03914': 'MNXM165365', 'C04427': 'MNXM169368', 'C04451': 'MNXM90601', 'C04721': 'MNXM94248', 'C05241': 'MNXM722784', 'C05997': 'MNXM36492', 'C07433': 'MNXM159617', 'C07970': 'MNXM58486', 'C08228': 'MNXM157170', 'C10528': 'MNXM1455', 'C11294': 'MNXM126764', 'C11458': 'MNXM724515', 'C11459': 'MNXM16413', 'C12995': 'MNXM527885', 'C13594': 'MNXM152253', 'C13595': 'MNXM157984', 'C14122': 'MNXM95692', 'C15870': 'MNXM1338', 'C18227': 'MNXM4582', 'C19718': 'MNXM1252', 'C80001': 'MNXM92466', 'C80002': 'MNXM4523', 'C80003': 'MNXM13824', 'C80004': 'MNXM13824', 'C80005': 'MNXM4219', 'C80006': 'MNXM13819', 'C80007': 'MNXM4214', 'C80008': 'MNXM32237', 'C80009': 'MNXM9521', 'C80010': 'MNXM9521', 'C80011': 'MNXM13821', 'C80012': 'MNXM73137', 'C80013': 'MNXM56307', 'C80014': 'MNXM4328', 'C80015': 'MNXM9763', 'C80016': 'MNXM14140', 'C80017': 'MNXM34598', 'C80018': 'MNXM14472', 'C80019': 'MNXM723734', 'C80020': 'MNXM34393', 'C80022': 'MNXM37357', 'C80025': 'MNXM14271', 'C80028': 'MNXM37315', 'C80029': 'MNXM14854', 'C80030': 'MNXM1801', 
                'C80031': 'MNXM6682', 'C80032': 'MNXM164200', 'C80034': 'MNXM2600', 'C80035': 'MNXM36033', 'C80036': 'MNXM6682', 'C80037': 'MNXM8046', 'C80038': 'MNXM107214', 'C80039': 'MNXM468570', 'C80040': 'MNXM468572', 'C80041': 'MNXM3800', 'C80042': 'MNXM2853', 'C80043': 'MNXM5084', 'C80045': 'MNXM11523', 'C80046': 'MNXM32005', 'C80047': 'MNXM78182', 'C80048': 'MNXM165478', 'C80049': 'MNXM165425', 'C80050': 'MNXM44968', 'C80051': 'MNXM165476', 'C80052': 'MNXM46488', 'C80053': 'MNXM73202', 'C80056': 'MNXM39260', 'C80057': 'MNXM14510', 'C80058': 'MNXM491041', 'C80061': 'MNXM5189', 'C80066': 'MNXM275', 'C80067': 'MNXM88927', 'C90010': 'MNXM31306', 'C90029': 'MNXM31741', 'C90040': 'MNXM31741', 'C90041': 'MNXM5401', 'C90073': 'MNXM1099', 'C90091': 'MNXM35167', 'C90092': 'MNXM164143', 'C90094': 'MNXM5476', 'C90096': 'MNXM35277', 'C90126': 'MNXM36486', 'C90131': 'MNXM162568', 'C90132': 'MNXM146349', 'C90134': 'MNXM36558', 'C90138': 'MNXM36756', 'C90139': 'MNXM36758', 'C90140': 'MNXM36762', 'C90144': 'MNXM36773', 'C90151': 'MNXM484827', 'C90162': 'MNXM37435', 'C90181': 'MNXM487477', 'C90205': 'MNXM38670', 'C90273': 'MNXM46158', 'C90316': 'MNXM411', 'C90323': 'MNXM411', 'C90325': 'MNXM507621', 'C90335': 'MNXM51449', 'C90336': 'MNXM11433', 'C90347': 'MNXM52367', 'C90348': 'MNXM52369', 'C90350': 'MNXM5738', 'C90352': 'MNXM52465', 'C90354': 'MNXM52466', 'C90366': 'MNXM53170', 'C90388': 'MNXM747', 'C90389': 'MNXM2377', 'C90390': 'MNXM1111', 'C90391': 'MNXM1112', 'C90411': 'MNXM56060', 'C90412': 'MNXM56121', 'C90413': 'MNXM724800', 'C90414': 'MNXM1190', 'C90418': 'MNXM56347', 'C90419': 'MNXM1067', 'C90427': 'MNXM57076 ', 'C90437': 'MNXM57486', 'C90438': 'MNXM527415', 'C90439': 'MNXM105743', 'C90444': 'MNXM57835', 'C90445': 'MNXM5189', 'C90497': 'MNXM59726', 'C90498': 'MNXM59727', 'C90534': 'MNXM63081', 'C90535': 'MNXM63174', 'C90536': 'MNXM63175', 'C90537': 'MNXM63157', 'C90538': 'MNXM63156', 'C90548': 'MNXM63755', 'C90549': 'MNXM63797', 'C90550': 'MNXM63798', 'C90552': 'MNXM64019', 'C90553': 'MNXM64020', 'C90563': 'MNXM64863', 'C90565': 'MNXM968', 'C90575': 'MNXM65303', 'C90579': 'MNXM811', 'C90590': 'MNXM66180', 'C90598': 'MNXM73155', 'C90599': 'MNXM7283', 'C90616': 'MNXM3337', 'C90640': 'MNXM81220', 'C90678': 'MNXM83846', 'C90694': 'MNXM114303', 'C90697': 'MNXM725810', 'C90703': 'MNXM5319 ', 'C90724': 'MNXM89469', 'C90725': 'MNXM89470', 'C90726': 'MNXM89471', 'C90727': 'MNXM2876', 'C90728': 'MNXM89472'}


    def _checkFilePath(self, path, filename):
        """Check that the directory and the filename are valid and choose to use
        either the local or the global path
        """
        if path and path[-1:]=='/':
            path = path[:-1]
        if path==None:
            if self.globalPath==None:
                logging.error('Both global path and local are not set')
                return None
            else:
                if os.path.isdir(self.globalPath):
                    if os.path.isfile(self.globalPath+'/'+filename):
                        return self.globalPath+'/'+filename
                    else:
                        logging.error('Global path file: '+str(filename)+', does not exist')
                        return None
                else:
                    logging.error('Global path is not a directory: '+str(self.globalPath))
                    return None
        else:
            if os.path.isdir(path):
                if os.path.isfile(path+'/'+filename):
                    return path+'/'+filename
                else:
                    logging.error('The file is not valid: '+str(path+'/'+filename))
                    return None
            else:
                logging.error('Local path is not a directory: '+str(path))
                return None


    def all(self):
        """Run all the files required to generate the cache
        Requires: mvc.db, chem_xref.tsv, chem_prop.tsv, cc_compounds.json.gz and alberty.json
        """
        if self.globalPath==None:
            logging.error('Need to define the global path to use all()')
            return False
        if not os.path.isdir(os.getcwd()+'/cache'):
            os.mkdir('cache')
        #deprecatedMNXM_mnxm
        logging.info('Generating deprecatedMNXM_mnxm')
        deprecatedMNXM_mnxm = self.deprecatedMNXM_mnxm()
        pickle.dump(deprecatedMNXM_mnxm, open('cache/deprecatedMNXM_mnxm.pickle', 'wb'))
        #mnxm_dG
        logging.info('Generating mnxm_dG')
        pickle.dump(self.mnxm_dG(deprecatedMNXM_mnxm), open('cache/mnxm_dG.pickle', 'wb'))
        #rr_reactions
        logging.info('Generating rr_reactions')
        pickle.dump(self.rr_reactions(), open('cache/rr_reactions.pickle', 'wb'))


    #TODO: save the deprecatedMNXM_mnxm to be used in case there rp_paths uses an old version of MNX
    def deprecatedMNXM_mnxm(self, chem_xref_path=None):
        """generate the deprecatedMNXM_mnxmdictionnary parameter from chem_xref.tsv from MNX
        """
        a = {}
        with open(self._checkFilePath(chem_xref_path, 'chem_xref.tsv')) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                mnx = row[0].split(':')
                if mnx[0]=='deprecated':
                    #a[row[1]] = mnx[1]
                    a[mnx[1]] = row[1]
            a.update(self.convertMNXM)
            a['MNXM01'] = 'MNXM1'
        return a


    def mnxm_dG(self,
                deprecatedMNXM_mnxm,
                chem_xref_path=None,
                cc_compounds_path=None,
                alberty_path=None):
        cc_alberty = {}
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
                            mnxm = deprecatedMNXM_mnxm[mnx[1]]
                        except KeyError:
                            mnxm = mnx[1]
                        kegg_mnxm[mnxm] = row[1]
        ####################### cc_compounds ############
        notFound_cc = []
        gz_file = gzip.open(self._checkFilePath(cc_compounds_path, 'cc_compounds.json.gz'), 'rb')
        f_c = gz_file.read()
        c = json.loads(f_c)
        for cd in c:
            try:
                mnxm = kegg_mnxm[cd['CID']]
            except KeyError:
                try:
                    mnxm = self.curated_kegg_mnxm[cd['CID']]
                    try:
                        mnxm = deprecatedMNXM_mnxm[mnxm]
                    except KeyError:
                        pass
                except KeyError:
                    logging.warning('Cannot find: '+str(cd))
                    notFound_cc.append(cd['CID'])
                    continue
            if not mnxm in cc_alberty:
                cc_alberty[mnxm] = {}
            if not 'component_contribution' in cc_alberty[mnxm]:
                cc_alberty[mnxm]['component_contribution'] = [cd]
            else:
                cc_alberty[mnxm]['component_contribution'].append(cd)
        ######################## alberty ################
        notFound_alberty = []
        with open(self._checkFilePath(alberty_path, 'alberty.json')) as json_data:
            d = json.loads(json_data.read())
        for cd in d:
            try:
                mnxm = kegg_mnxm[cd['cid']]
            except KeyError:
                try:
                    mnxm = self.curated_kegg_mnxm[cd['cid']]
                    try:
                        mnxm = deprecatedMNXM_mnxm[mnxm]
                    except KeyError:
                        pass
                except KeyError:
                    logging.warning('Cannot find: '+str(cd))
                    notFound_alberty.append(cd['cid'])
                    continue
            if not mnxm in cc_alberty:
                cc_alberty[mnxm] = {}
            if not 'alberty' in cc_alberty[mnxm]:
                cc_alberty[mnxm]['alberty'] = [cd]
            else:
                cc_alberty[mnxm]['alberty'].append(cd)
        return cc_alberty


    #implement parsing a MNX flat file for the generation of the rr_reactions flat file
    #TODO: change the name of the parameter to generate the correct one
    def reactions_mnx(self):
        return False


    #DO NOT DELETE -- to generate rr_reactions.pickle
    #WARNING: there could be instances where two rules are replaced (CHECK)
    def rr_reactions(self, mvc_path=None):
        """Generate the rr_reactions.pickle from the RetroRules database
        Function that parses retrorules to generate rp2paths ID with the metanetX ID's for the reactions
        Used to be pickled and opened upon initiation of the parser object
        """
        ################# extract all the products from RetroRules ####################
        try:
            ########## open either the global path or the local defined path ############
            #### (with priority with the local path)
            conn = None
            db_file = self._checkFilePath(mvc_path, 'mvc.db')
            if db_file:
                conn = sqlite3.connect(db_file)
            else:
                return {}
            c = conn.cursor()
            ################### select all the unique rules ID ########
            c.execute('''SELECT DISTINCT reactions.mnxr, chemical_species.mnxm FROM rules
                        INNER JOIN reactions 
                            ON (reactions.id=rules.reaction_id)
                        INNER JOIN chemical_species 
                            ON (chemical_species.id=rules.substrate_id) ''')
            all_rows = c.fetchall()
            reactions = {}
            rule_ids = {}
            for row in all_rows:
                reactions[row[0]+'_'+row[1]] = {'right': {}, 'left': {}} #unique
                if not row[0] in rule_ids: #possibly not unique
                    rule_ids[row[0]] = []
                #not sure about this TODO
                mnxm = row[1]
                if mnxm in self.convertMNXM:
                    mnxm = self.convertMNXM[mnxm]
                rule_ids[row[0]].append(row[0]+'_'+mnxm)
            ################### select products for each given reaction #############
            c.execute(''' SELECT reactions.mnxr,
                                chemical_species.mnxm,
                                reaction_products.stochiometry 
                                    FROM reaction_products
                            INNER JOIN reactions 
                                ON (reactions.id=reaction_products.reaction_id)
                            INNER JOIN chemical_species
                                ON (chemical_species.id=reaction_products.chemical_id)  ''')
            all_rows = c.fetchall()
            for row in all_rows:
                if not row[0] in rule_ids:
                    rule_ids[row[0]] = {}
                for rule in rule_ids[row[0]]:
                    toAdd = row[1]
                    if toAdd in self.convertMNXM:
                        toAdd = self.convertMNXM[toAdd]
                    reactions[rule]['right'][toAdd] = int(row[2])
            ################### select substrates for each given reaction #############
            c.execute(''' SELECT reactions.mnxr,
                                chemical_species.mnxm,
                                reaction_substrates.stochiometry 
                                    FROM reaction_substrates
                            INNER JOIN reactions 
                                ON (reactions.id=reaction_substrates.reaction_id)
                            INNER JOIN chemical_species
                                ON (chemical_species.id=reaction_substrates.chemical_id)  ''')
            all_rows = c.fetchall()
            for row in all_rows:
                if not row[0] in rule_ids:
                    rule_ids[row[0]] = {}
                for rule in rule_ids[row[0]]:
                    toAdd = row[1]
                    if toAdd in self.convertMNXM:
                        toAdd = self.convertMNXM[toAdd]
                    reactions[rule]['left'][toAdd] = int(row[2])
            #self.rr_reactions = reactions
            return reactions 
        except sqlite3.OperationalError:
            logging.error('Error parsing the sqlite3 database ('+str(path)+')')
            conn.close()
            return {}

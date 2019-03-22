import csv
import os
import itertools
import cobra
import logging
import collections
import numpy as np
import pickle
import copy
import gzip
import libsbml
import sqlite3
import re
from cameo.strain_design import OptKnock
from hashlib import md5
from cobra import Model, Reaction, Metabolite



class SBML:
    """docstring for ClassName"""
    def __init__(self, outputPath=None, Id_type=None):
        '''if inputPath and inputPath[-1:]=='/':
            inputPath = inputPath[:-1]'''
        if outputPath and outputPath[-1:]=='/':
            outputPath = outputPath[:-1]
        self.mvc = None
        self.outputPath = outputPath
        self.Id_type = Id_type

        


    ########################### PRIVATE FUNCTIONS ###############################


    #TODO: export the models into SBML, controlling the versions
    #TODO: export the heterologous reactions into SBML format
    def _exportSBML(self, type_name, model, model_id, inputType, path=None):
        """
        Function to export models generated in libSBML or cobraPy to SBML
        cobra.io.write_legacy_sbml(models[model_id], 
                                    p+'/'+str(model_id)+'.sbml', 
                                    use_fbc_package=False)
        """
        ####### check the path #########
        p = None
        if path:
            if path[-1:]=='/':
                path = path[:-1]
            if not os.path.isdir(path):
                if self.outputPath:
                    p = self.outputPath
                else:
                    logging.error('The output path is not a directory: '+str(path))
                    return False
            else:
                p = path
        else:
            p = self.outputPath
        p += '/'+type_name
        ########## check and create folder #####
        if not os.path.exists(p):
            os.makedirs(p)
        ########## export ###################
        """
        for model_id in models:
            if inputType=='cobrapy':
                cobra.io.write_sbml_model(models[model_id],
                                            p+'/'+str(model_id)+'.sbml', 
                                            use_fbc_package=False)
            elif inputType=='libsbml':
                libsbml.writeSBMLToFile(models[model_id],
                                p+'/'+str(model_id)+'.sbml')
            else:
                logging.error('Cannot recognise input type')
                return False
        """
        if inputType=='cobrapy':
            cobra.io.write_sbml_model(model,
                                        p+'/'+str(model_id)+'.sbml', 
                                        use_fbc_package=True)
        elif inputType=='libsbml':
            libsbml.writeSBMLToFile(model,
                            p+'/'+str(model_id)+'.sbml')
        else:
            logging.error('Cannot recognise input type')
            return False
        return True
    ########################### PUBLIC FUNCTIONS ###############################

    def all(self, cofactors_rp_paths, dico_mvc, model_compartments, cobra_model, path=None, isCPLEX=False, isExport=False):
        self.mvc = dico_mvc
        for path_id in cofactors_rp_paths:
            logging.info('################## SBML creation for path '+str(path_id)+'  #####################')
            tmp_model = self.constructModel_cobraPy(cofactors_rp_paths[path_id], 
                                                    path_id, 
                                                    model_compartments, 
                                                    cobra_model, 
                                                    isExport=True)
            if isCPLEX:
                tmp_model.solver = 'cplex'
        return True
    

    #Given the path, open the model (NOTE: only MNX models for now)
    def model(self, path=None):
        try:
            return cobra.io.read_sbml_model(self._checkFilePath(path, 'model'))
        except AttributeError:
            return None

    #Function to change the id of the metabolite to metanetx and the user's one
    def change_id(self, meta, number_id, mvc_path="input_cache/"):
        liste = self.mvc
        check = True
        for i in liste:
            for x,z in enumerate(i.split(':')):
                if meta == z:
                    if i.split(':')[number_id] == 'None':
                        logging.info('Identifier not found for compound '+meta+', provided MetaNetX one.')
                        tmp_id = i.split(':')[0]
                    else:
                        tmp_id = i.split(':')[number_id]
                    check = False
            if not check:
                break
        return tmp_id

    #def constructModel_cobrapy(self, model_id, model_compartments, isExport=False, inPath=None):
    def constructModel_cobraPy(self, cofactors_rp_path, model_id, ori_model_compartments, ori_model, isExport=False, inPath=None):
        
        #Defining the identifier to be choosen in the list
        Id = self.Id_type
        if Id in ('metanetx', 'MetaNetx', 'MetanetX'):
            Id_ref = 0
        elif Id in ('kegg', 'Kegg'):
            Id_ref = 1
        elif Id in ('seed', 'Seed'):
            Id_ref = 2
        elif Id in ('metacyc', 'Metacyc', 'MetaCyc'):
            Id_ref = 3
        elif Id in ('reactome', 'Reactome'):
            Id_ref = 4
        elif Id in ('sabiork', 'Sabiork'):
            Id_ref = 5
        elif Id in ('bigg', 'Bigg'):
            Id_ref = 6
        elif Id in ('Chebi', 'chebi'):
            Id_ref = 7
        elif Id in ('Hmdb', 'hmdb', 'HMDB'):
            Id_ref = 8

        #smiles = self.rp_smiles
        
        """
            Returns a dictionnary of models with the keys the path from RP2paths out_paths.csv
            and the instructions to plot the heterologous pathway (as well as the metabolic sink and the source)
            to generate a plot in networkx
        """
        cytoplasm_compartment = [i for i in ori_model_compartments if ori_model_compartments[i]['short_name']=='c'][0] # this assuming that there is always only one result
        extracellular_compartment = [i for i in ori_model_compartments if ori_model_compartments[i]['short_name']=='e'][0] # this assuming that there is always only one result
        #list all the metabolites in the model. WARNING: works only for MNXM models
        #all_inputModel_metabo = [i.id.split('__')[0] for i in ori_model.metabolites]
        #NOTE: we assume that the metabolites are all in the cytoplasm (apart from the last transport step)
        #TODO: need to flag that the metabolites (sink) first step in the reaction is contained in the model - to validate the rp_path
        ########### METABOLITES #########################
        #create a new model where we will add this path to it
        
        
        #enumerate all the different compounds from the path
        #new_meta = list(set([y for i in path for y in itertools.chain(i['right'], i['left'])]))
        new_meta = set([i for step_id in cofactors_rp_path['path'] for i in cofactors_rp_path['path'][step_id]['steps'].keys()])
        
        if ori_model == None:
            model = Model('%s' % model_id)
        else :
            model = ori_model.copy()

        def avecModel():
            #remove the ones that already exist in the model
            for meta in new_meta:
                if not meta in all_meta:
                    try:
                        #NOTE: we assume that all the metabolites that we add here are in the cytoplasm
                        #return the metaolite if already in the model
                        all_meta[meta] = model.metabolites.get_by_id(meta+'__64__'+cytoplasm_compartment)
                        tmp_meta = self.change_id(meta, Id_ref) ##!!!!!!
                        #tmp_meta = model.metabolites.get_by_id(meta.split(':')[0]+'__64__'+cytoplasm_compartment)
                        all_meta[meta].id = "%s" %  tmp_meta
                        all_meta[meta].compartment = "c"
                        #all_meta[meta.split(':')[0]] = tmp_meta
                    except KeyError:
                        #if not in the model create a new one
                        logging.info('Compound not in the model, creating model metabolite for '+meta.split(':')[0])
                        all_meta[meta] = cobra.Metabolite(meta, name=meta, compartment=cytoplasm_compartment)
                        if 'MNXM' in meta[:4]:    
                            tmp_meta = self.change_id(meta, Id_ref) ##!!!!!!
                            all_meta[meta].id = "%s" %  tmp_meta
                            all_meta[meta].compartment = "c"
            return all_meta[meta]

        def sansModel():
            for meta in new_meta:
                if not meta in all_meta:
                    all_meta[meta] = cobra.Metabolite(meta, name=meta, compartment='cytoplasm')
                    if "MNXM" in meta[:4]:  # [TODO] change it to if not [CMPD] because i will not use the MNXM identifiers
                        tmp_meta = self.change_id(meta, Id_ref)
                        all_meta[meta].id = "%s" %  tmp_meta
            return all_meta[meta]

        switcher = {
                ori_model: avecModel,
                None: sansModel
                }

        all_meta = {}
        
        func = switcher.get(ori_model, "Error choosing option")
        func()
       

        ############## REACTIONS ##########################
        for step_id in cofactors_rp_path['path']:
            reaction = cobra.Reaction('rpReaction_'+str(step_id))
            #reaction.name = 
            reaction.lower_bound = 0.0 # assume that all the reactions are irreversible
            reaction.upper_bound = 999999.0 #this is dependent on the fluxes of the others reactions
            reaction.gene_reaction_rule = 'rpGene_'+str(step_id)
            reac_meta = {}
            for mnxm in cofactors_rp_path['path'][step_id]['steps']:
                reac_meta[all_meta[mnxm]] = float(cofactors_rp_path['path'][step_id]['steps'][mnxm]['stoichiometry']) #mnxm.split(':')[0]
            reaction.add_metabolites(reac_meta)
            model.add_reactions([reaction]) #just model
        ################# Extracellular transport of target
        #NOTE: some molecules cannot be exported and thus this step should not be added
        #identify the target molecule
        target_name = [i for i in new_meta if i[:6]=='TARGET'][0]
        #create the metabolite for the extracellular version of the target metabolite
        extracell_target = cobra.Metabolite(target_name+'_e', name=target_name, compartment=extracellular_compartment)
        #Add the export from the cytoplasm to the extracellular matrix ######
        #exportReaction = cobra.Reaction('exportTarget')
        #exportReaction.name = 'ExportTarget'
        exportReaction = cobra.Reaction('targetSink')
        exportReaction.name = 'targetSink'
        exportReaction.lower_bound = 0.0 #default = 0.0
        exportReaction.upper_bound = 999999.0 #default = 1000 TODO: see if changing this changes something
        #these are the bounds for the yeast bigg model
        #TODO: check .reversibility of the reaction after setting these bounds
        #add that metabolite to the 
        #exportReaction.add_metabolites({extracell_target: 1.0, all_meta[target_name]: -1.0})
        exportReaction.add_metabolites({all_meta[target_name]: -1.0}) # target_name.split(':')[0]
        model.add_reactions([exportReaction])
        ################### Add the sink reaction
        '''
        sinkReaction = cobra.Reaction('targetSink')
        sinkReaction.name = 'TargetSink'
        sinkReaction.lower_bound = 0.0 # we assume that all the reactions are irreversible
        sinkReaction.upper_bound = 999999.0 # this is dependent on the fluxes of the other reactions
        sinkReaction.add_metabolites(
            {extracell_target: -1.0})
        model.add_reactions([sinkReaction])
        '''
        if isExport:
            self._exportSBML('sbml_models', model, model_id, 'cobrapy', inPath)
        return model
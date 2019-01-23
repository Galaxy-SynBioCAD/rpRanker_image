import cobra
import itertools
import os
from cameo.strain_design import OptKnock

class FBA:
    def __init__(self, outputPath=None):
        if outputPath and outputPath[-1:]=='/':
            outputPath = outputPath[:-1]
        self.outputPath = outputPath
        self.rp_sbml_paths = None
        self.rp_sbml_models = None
        self.results = None


    ############################### PRIVATE FUNCTIONS ############################


    def _switchToCPLEX(self, models):
        """ 
            All the models are swiched to using the CPLEX solver
        """
        #TODO: have a check to see if you can access the CPLEX solver before assigning it to each model
        for model_id in models:
            models[model_id].solver = 'cplex'


    #TODO: export the models into SBML, controlling the versions
    #TODO: export the heterologous reactions into SBML format
    def _exportSBML(self, type_name, models, path=None):
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
        for model_id in models:
            cobra.io.write_sbml_model(models[model_id], 
                                        p+'/'+str(model_id)+'.sbml', 
                                        use_fbc_package=False)
            '''
            cobra.io.write_legacy_sbml(models[model_id], 
                                        p+'/'+str(model_id)+'.sbml', 
                                        use_fbc_package=False)
            '''
        return True


    def _check(value, message):
        """
        Taken from: http://sbml.org/Software/libSBML/docs/python-api/create_simple_model_8py-example.html
        """
        if value is None:
            raise SystemExit('LibSBML returned a null value trying to ' + message + '.')
        elif type(value) is int:
            if value == LIBSBML_OPERATION_SUCCESS:
                return
            else:
                err_msg = 'Error encountered trying to ' + message + '.' \
                        + 'LibSBML returned error code ' + str(value) + ': "' \
                        + OperationReturnValue_toString(value).strip() + '"'
                raise SystemExit(err_msg)
        else:
            return


    ############################### PUBLIC FUNCTIONS ############################# 

    #function that takes the location of the database, the model of interest, and the RetroPath2.0 output
    #and constructs models including the heterologous pathway
    def runAll(self, cofactors_rp_paths, rr_reactions, model_compartments, cobra_model, rp_smiles, isCPLEX=False):
        """run all the FBA models
            Adds the cofactors to the heterologous pathways and then to the models 
            and runs the models using different objectives
        """
        self.rp_sbml_paths = self.constructPaths(cofactors_rp_paths, model_compartments, rp_smiles, True)
        self.rp_sbml_models = self.constructModels(cofactors_rp_paths, model_compartments, cobra_model, True)
        if isCPLEX:
            self._switchToCPLEX(self.rp_sbml_models)
        self.results = {}
        '''
        self.results['biomass'] =  self.simulateBiomass(self.rp_sbml_models)
        self.results['target'] = self.simulateTarget(self.rp_sbml_models)
        self.results['biLevel'] = self.simulateBiLevel(self.rp_sbml_models)
        self.results['splitObjective'] = self.simulateSplitObjective(self.rp_sbml_models)
        '''


        self.results['biomass'] = {}
        self.results['biomass']['cameo'], self.results['biomass']['sorted'] = self.simulateBiomass(self.rp_sbml_models)
        self.results['target'] = {}
        self.results['target']['cameo'], self.results['target']['sorted'] = self.simulateTarget(self.rp_sbml_models)
        self.results['splitObjective'] = {}
        self.results['splitObjective']['cameo'], self.results['splitObjective']['sorted'] = self.simulateSplitObjective(self.rp_sbml_models)
        self.results['biLevel'] = {}
        self.results['biLevel']['cameo'], self.results['biLevel']['sorted'] = self.simulateBiLevel(self.rp_sbml_models)
        return True 


    def constructModels(self, cofactors_rp_paths, ori_model_compartments, ori_model, isExport=False, inPath=None):
        """
            Returns a dictionnary of models with the keys the path from RP2paths out_paths.csv
            and the instructions to plot the heterologous pathway (as well as the metabolic sink and the source)
            to generate a plot in networkx
        """
        all_rp_models = {}
        cytoplasm_compartment = [i for i in ori_model_compartments if ori_model_compartments[i]['short_name']=='c'][0] # this assuming that there is always only one result
        extracellular_compartment = [i for i in ori_model_compartments if ori_model_compartments[i]['short_name']=='e'][0] # this assuming that there is always only one result
        #list all the metabolites in the model. WARNING: works only for MNXM models
        all_inputModel_metabo = [i.id.split('__')[0] for i in ori_model.metabolites]
        #NOTE: we assume that the metabolites are all in the cytoplasm (apart from the last transport step)
        #TODO: need to flag that the metabolites (sink) first step in the reaction is contained in the model - to validate the rp_path
        for path_id in cofactors_rp_paths:
            ########### METABOLITES #########################
            #create a new model where we will add this path to it
            model = ori_model.copy()
            #enumerate all the different compounds from the path
            #new_meta = list(set([y for i in path for y in itertools.chain(i['right'], i['left'])]))
            new_meta = set([i for path_id in cofactors_rp_paths for step_id in cofactors_rp_paths[path_id]['path'] for i in cofactors_rp_paths[path_id]['path'][step_id]['step'].keys()])
            all_meta = {}
            for meta in new_meta:
                #remove the ones that already exist in the model
                if not meta in all_meta:
                    try:
                        #NOTE: we assume that all the metabolites that we add here are in the cytoplasm
                        #return the metaolite if already in the model
                        all_meta[meta] = model.metabolites.get_by_id(meta+'__64__'+cytoplasm_compartment)
                    except KeyError:
                        #if not in the model create a new one
                        all_meta[meta] = cobra.Metabolite(meta, name=meta, compartment=cytoplasm_compartment)
            ############## REACTIONS ##########################
            for step_id in cofactors_rp_paths[path_id]['path']:
                reaction = cobra.Reaction('rpReaction_'+str(step_id))
                #reaction.name = 
                reaction.lower_bound = 0.0 # assume that all the reactions are irreversible
                reaction.upper_bound = 999999.0 #this is dependent on the fluxes of the others reactions
                reaction.gene_reaction_rule = 'rpGene_'+str(step_id)
                reac_meta = {}
                for mnxm in cofactors_rp_paths[path_id]['path'][step_id]['step']:
                    reac_meta[all_meta[mnxm]] = float(cofactors_rp_paths[path_id]['path'][step_id]['step'][mnxm]['stochio'])
                reaction.add_metabolites(reac_meta)
                model.add_reactions([reaction])
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
            exportReaction.add_metabolites({all_meta[target_name]: -1.0})
            #print(exportReaction.data_frame())  
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
            all_rp_models[path_id] = model
        if isExport:
            self._exportSBML('sbml_models', all_rp_models, inPath)
        return all_rp_models




    def constructPaths_libSBML():
        """
        Function to construct a series of libSBML objects for using FBC (constraint based) 
        """
        for path_id in cofactors_rp_paths:
            sbmlDoc = libsbml.SBMLDocument(3,2) #level, version
            #Model
            model = libsbml.createModel()
            self._check(model, 'create model')
            self._check(model.setTimeUnits('second'), 'set model-wide time units')
            self._check(model.setExtentUnits('mole'), 'set model units of extent')
            self._check(model.setSubstanceUnits('mole'), 'set model substance units')
            model.setId('Path '+str(path_id))
            #Unit definition
            per_second = model.createUnitDefinition()
            self._check(per_second, 'create unit definition')
            self._check(per_second.setId('per_second'), 'set unit definition id')
            unit = per_second.createUnit()
            self._check(unit, 'create unit on per_second')
            self._check(unit.setKind(UNIT_KIND_SECOND), 'set unit kind')
            self._check(unit.setExponent(-1), 'set unit exponent')
            self._check(unit.setScale(0), 'set unit scale')
            self._check(unit.setMultiplier(1), 'set unit multiplier')
            #Compartments
            cytoplasm = model.createCompartment()
            self._check(cytoplasm, 'create compartment')
            self._check(cytoplasm.setId('cytoplasm'), 'set compartment id')
            self._check(cytoplasm.setConstant(True), 'set compartment "constant"')
            self._check(cytoplasm.setSize(1), 'set compartment "size"')
            self._check(cytoplasm.setSpatialDimensions(3), 'set compartment dimensions')
            self._check(cytoplasm.setUnits('litre'), 'set compartment size units')
            #Species
            new_meta = set([i for path_id in cofactors_rp_paths for step_id in cofactors_rp_paths[path_id]['path'] for i in cofactors_rp_paths[path_id]['path'][step_id]['step'].keys()])
            for meta in new_meta:
                #TODO: annotate it with MNXC3 id for cytoplasm compartment
                #Species --> Loop through all the RP paths
                spe = model.createSpecies()
                self._check(spe, 'create species spe')
                self._check(spe.setId('spe'), 'set species spe id')
                self._check(spe.setCompartment('c1'), 'set species spe compartment')
                self._check(spe.setConstant(False), 'set "constant" attribute on spe')
                self._check(spe.setInitialAmount(1), 'set initial amount for spe')
                self._check(spe.setSubstanceUnits('mole'), 'set substance units for spe')
                self._check(spe.setBoundaryCondition(False), 'set "boundaryCondition" on spe')
                self._check(spe.setHasOnlySubstanceUnits(False), 'set "hasOnlySubstanceUnits" on spe')
            #Reactions
            for step_id in cofactors_rp_paths[path_id]['path']:
                r1 = model.createReaction()
                check(r1, 'create reaction')
                check(r1.setId('r1'), 'set reaction id')
                check(r1.setReversible(False), 'set reaction reversibility flag')
                check(r1.setFast(False), 'set reaction "fast" attribute')
                species_ref1 = r1.createReactant()
                check(species_ref1, 'create reactant')
                check(species_ref1.setSpecies('s1'), 'assign reactant species')
                check(species_ref1.setConstant(True), 'set "constant" on species ref 1')
                species_ref2 = r1.createProduct()
                check(species_ref2, 'create product')
                check(species_ref2.setSpecies('s2'), 'assign product species')
                check(species_ref2.setConstant(True), 'set "constant" on species ref 2')
                #Annotations for reaction
                cv = libsbml.CVTerm()
                cv.setQualifierType(libsbml.BIOLOGICAL_QUALIFIER)
                cv.setBiologicalQualifierType(libsbml.BQB_IS_DESCRIBED_BY)
                cv.addResource("https://retrorules.org/RULE/RR00121121/16/False")
                react.addCVTerm(cv)





    def constructPaths(self, cofactors_rp_paths, ori_model_compartments, rp_smiles, isExport=False, inPath=None):
        """Construct the cobra models of the RetroPath paths with the cofactors and export them to SBML

            Assumes that all the reactions are happening in the cytoplasm
        """
        rp_sbml_paths = {}
        cytoplasm_compartment = [i for i in ori_model_compartments if ori_model_compartments[i]['short_name']=='c'][0] # this assuming that there is always only one result
        extracellular_compartment = [i for i in ori_model_compartments if ori_model_compartments[i]['short_name']=='e'][0] # this assuming that there is always only one result
        #NOTE: we assume that the metabolites are all in the cytoplasm (apart from the last transport step)
        #TODO: need to flag that the metabolites (sink) first step in the reaction is contained in the model - to validate the rp_path
        for path_id in cofactors_rp_paths:
            ########### METABOLITES #########################
            #create a new model where we will add this path to it
            model = cobra.Model(str(path_id))
            #new_meta=list(set([y for i in cofactors_rp_paths[path] for y in itertools.chain(i['right'], i['left'])]))
            new_meta = set([i for path_id in cofactors_rp_paths for step_id in cofactors_rp_paths[path_id]['path'] for i in cofactors_rp_paths[path_id]['path'][step_id]['step'].keys()])
            all_meta = {}
            #enumerate all the unique compounds from a path
            for meta in new_meta:
                #remove the ones that already exist in the model
                if not meta in all_meta:
                    all_meta[meta] = cobra.Metabolite(meta, name=meta, compartment=cytoplasm_compartment)
                    all_meta[meta].annotation = { 'ImaginaryCompDB':'SpecificCompIdentifier', "uniprot": ["Q12345", "P12345"]}
            ############## REACTIONS ##########################
            for step_id in cofactors_rp_paths[path_id]['path']:
                #TODO: could this lead to a conflict if there is the same reactionID and a different substrateID
                #reaction = cobra.Reaction(step['rule_id'].split('_')[0])
                reaction = cobra.Reaction('rpReaction_'+str(step_id))
                reaction.lower_bound = 0.0 # assume that all the reactions are irreversible
                reaction.upper_bound = 999999.0 #this is dependent on the fluxes of the others reactions
                reaction.gene_reaction_rule = 'RPGene_'+str(step_id)
                reac_meta = {}
                for mnxm in cofactors_rp_paths[path_id]['path'][step_id]['step']:
                    reac_meta[all_meta[mnxm]] = float(cofactors_rp_paths[path_id]['path'][step_id]['step'][mnxm]['stochio'])
                reaction.add_metabolites(reac_meta)
                model.add_reactions([reaction])
            ################# Extracellular transport of target
            #NOTE: some molecules cannot be exported and thus this step should not be added
            #identify the target molecule
            target_name = [i for i in new_meta if i[:6]=='TARGET'][0]
            #create the metabolite for the extracellular version of the target metabolite
            extracell_target = cobra.Metabolite(target_name+'_e',
                                                name = target_name,
                                                compartment = extracellular_compartment)
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
            exportReaction.add_metabolites({all_meta[target_name]: -1.0})
            #print(exportReaction.data_frame())  
            model.add_reactions([exportReaction])
            #rp_sbml_paths[path[0]['path_id']] = model
            rp_sbml_paths[path_id] = model
        if isExport:
            self._exportSBML('sbml_paths', rp_sbml_paths, inPath)
        return rp_sbml_paths


    ############################### FBA pathway ranking ####################

    #1) Number of interventions

    #2) Maximal growth rate

    def simulateBiomass(self, all_rp_models):
        biomass = {}
        for model_id in all_rp_models:
            biomass[model_id] = {}
            biomass[model_id] = all_rp_models[model_id].optimize()
        #return biomass, sorted(biomass.items(), key=lambda x: x[1], reverse=True)
        return biomass, sorted([(biomass[i].fluxes['targetSink'], i) for i in biomass.keys()], reverse=True)
        #TODO: return a sorted list of the best performing pathways


    #3) Minimum product yeild at maximal growth rate

    #4) Minimum product yeild

    #5) Anaerobic condition

    #6) Number of potentially disruptive products

    #7) Number of accessible metabolites (avoid intermediate accumulation)

    #8) Thermodynamics (MDF)

    #9) The overlap of the same changes --> might not be applicable in our case

    #10) Reduced model:q


    def simulateTarget(self, all_rp_models):
        target = {}
        for model_id in all_rp_models:
            target[model_id] = {}
            all_rp_models[model_id].objective = 'targetSink'
            target[model_id] = all_rp_models[model_id].optimize()
        #return target, sorted(target.items(), key=lambda x: x[1], reverse=True)
        return target, sorted([(target[i].fluxes['targetSink'], i) for i in target.keys()], reverse=True)
        #TODO: return a sorted list of the best performing pathways


    def simulateSplitObjective(self, all_rp_models, biomass=0.5):
        if not 0.0<biomass<1.0:
            logging.error('The proportion of the objective that is given to BIOMASS must be 0.0< and 1.0>')
            return {}
        splitObj = {}
        for model_id in all_rp_models:
            #TODO: need to define what is the biomass reaction from function input
            all_rp_models[model_id].objective = {all_rp_models[model_id].reactions.R48E37591: biomass, all_rp_models[model_id].reactions.targetSink: 1.0-biomass}
            splitObj[model_id] = all_rp_models[model_id].optimize()
        #return split_Obj, sorted(splitObj.items(), key=lambda x: x[1], reverse=True)
        return splitObj, sorted([(splitObj[i].fluxes['targetSink'], i) for i in splitObj.keys()], reverse=True)
        #TODO: return a sorted list of the best performing pathways


    def simulateBiLevel(self, all_rp_models):
        sim_res = {}
        for model_id in all_rp_models:
            sim_res[model_id] = {}
            try:
                #parameter exclude_non_gene_reactions --> 
                #   since we are introducing new pathways that are indeed not associated
                #   with genes we need to remove this parameter
                optknock = OptKnock(all_rp_models[model_id], exclude_non_gene_reactions=False, remove_blocked=False)
                sim_res[model_id] = optknock.run(max_knockouts=0, target='targetSink', biomass='R48E37591')
            except KeyError:
                print('KeyError with targetSink.... Not sure why that is')
        #return sim_res, sorted(sim_res.items(), key=lambda x: x[1], reverse=True)
        return sim_res, sorted([(sim_res[i].fluxes[0]['targetSink'], i) for i in sim_res.keys()], reverse=True)



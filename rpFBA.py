import cobra
import lxml
import rpSBML

class rpFBA:
    def __init__(self, rpsbml):
        self.rpsbml = rpsbml

    def libsbml_to_cobra(self):
        return cobra.io.sbml3.parse_xml_into_model(lxml.etree.fromstring(rpsbml.document.toXMLNode().toXMLString()))

    



import cobra
import itertools
import os
from cameo.strain_design import OptKnock
import libsbml
from hashlib import md5
import logging




#TODO: define the biomass parameter here to be used. The default should be the BIGG (R48E37591)
class FBA:
    def __init__(self, outputPath=None, biomassID='R48E37591'):
        if outputPath and outputPath[-1:]=='/':
            outputPath = outputPath[:-1]
        self.outputPath = outputPath
        self.biomassID = biomassID
        self.rp_sbml_paths = None
        self.rp_sbml_models = None
        self.results = None


    ########################################################################
    ############################### FBA pathway ranking ####################
    ########################################################################

    #1) Number of interventions
    # need to calculate the number of steps that are not native to know the number of interventions

    #2) Maximal growth rate

    def simulateBiomass(self, cofactors_rp_path, model):
        #TODO: update the objective function here
        res = model.optimize()
        cofactors_rp_path['flux_biomass'] = res.fluxes['targetSink']
        for step_id in cofactors_rp_path['path']:
            cofactors_rp_path['path'][step_id]['flux_biomass'] = res.fluxes['rpReaction_'+str(step_id)]
        return True

    #3) Minimum product yeild at maximal growth rate

    #4) Minimum product yeild

    #5) Anaerobic condition

    #6) Number of potentially disruptive products

        #Toxicity?

    #7) Number of accessible metabolites (avoid intermediate accumulation)

    #8) Thermodynamics (MDF)

    #9) The overlap of the same changes --> might not be applicable in our case

    #10) Reduced model

    #11) ECM

    def simulateTarget(self, cofactors_rp_path, model):
        model.objective = 'targetSink'
        res = model.optimize()
        cofactors_rp_path['flux_target'] = res.fluxes['targetSink']
        for step_id in cofactors_rp_path['path']:
            cofactors_rp_path['path'][step_id]['flux_target'] = res.fluxes['rpReaction_'+str(step_id)]
        return True


    def simulateSplitObjective(self, cofactors_rp_path, model, ratio_biomass=0.5):
        if not 0.0<ratio_biomass<1.0:
            logging.error('The proportion of the objective that is given to BIOMASS must be 0.0< and 1.0>')
            return False
        model.objective = {model.reactions.R48E37591: ratio_biomass, 
                            model.reactions.targetSink: 1.0-ratio_biomass}
        res = model.optimize()
        cofactors_rp_path['flux_splitObj'] = res.fluxes['targetSink']
        for step_id in cofactors_rp_path['path']:
            cofactors_rp_path['path'][step_id]['flux_splitObj'] = res.fluxes['rpReaction_'+str(step_id)] 
        return True


    def simulateBiLevel(self, cofactors_rp_path, model):
        try:
            optknock = OptKnock(model, exclude_non_gene_reactions=False, remove_blocked=False)
            sim_res = optknock.run(max_knockouts=0, target='targetSink', biomass=self.biomassID)
            cofactors_rp_path['flux_biLevel'] = sim_res.fluxes[0]['targetSink']
        except KeyError:
            logging.error('KeyError with targetSink.... Not sure why that is')
            return False
        for step_id in cofactors_rp_path['path']:
            cofactors_rp_path['path'][step_id]['flux_biLevel'] = sim_res.fluxes[0]['rpReaction_'+str(step_id)]
        return True


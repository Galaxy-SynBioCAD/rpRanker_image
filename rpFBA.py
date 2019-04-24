import cobra
import lxml
import rpSBML


## Class to simulate an rpsbml object using different FBA types and objective functions
#
# At this point we want to have the BIOMASS, target and shared objective
class rpFBA:
    def __init__(self, rpsbml):
        self.rpsbml = rpsbml
        self.cobraModel = None
  

    ## Given that there can be nultiple objectives defined in the SBML, this function switches between different ones
    #
    # TODO: 
    def switchObjective(self, objId):
        fbc_plugin = rpsbml.model.getPlugin('fbc')
        listObj = fbc_plugin.getListOfObjectives()
        if objId in [i.getId() for i in listObj]:
            listObj.setActiveObjective(objId)
        else:
            logger.warning('The objective Id '+str(objId)+' does not exist in rpsbml')
        self.libsbml_to_cobra()


    def libsbml_to_cobra(self):
        self.cobraModel = cobra.io.read_sbml_model(document.toXMLNode().toXMLString())
        #for an old version of cobrapy (0.4)
        #return cobra.io.sbml3.parse_xml_into_model(lxml.etree.fromstring(rpsbml.document.toXMLNode().toXMLString()))
        

    def simulate(self):
        res = self.cobraModel.optimize()
        #TODO: find a way of storing this information 



    ########################################################################
    ############################### FBA pathway ranking ####################
    ########################################################################

    #1) Number of interventions
    # need to calculate the number of steps that are not native to know the number of interventions

    #2) Maximal growth rate


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



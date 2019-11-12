import cobra
import libsbml

import logging
#import rpSBML

## Class to simulate an rpsbml object using different FBA types and objective functions
#
# At this point we want to have the BIOMASS, target and shared objective
#TODO: add the pareto frontier optimisation as an automatic way to calculate the optimal fluxes
class rpFBA:
    def __init__(self, rpsbml):
        self.logger = logging.getLogger(__name__)
        self.logger.info('Started instance of rpFBA')
        self.rpsbml = rpsbml
        self.cobraModel = None


    ## Check the libSBML calls
    #
    # Check that the libSBML python calls do not return error INT and if so, display the error. Taken from: http://sbml.org/Software/libSBML/docs/python-api/create_simple_model_8py-example.html
    #
    # @param value The SBML call
    # @param message The string that describes the call
    def _checklibSBML(self, value, message):
        if value is None:
            self.logger.error('LibSBML returned a null value trying to ' + message + '.')
            raise SystemExit('LibSBML returned a null value trying to ' + message + '.')
        elif type(value) is int:
            if value==libsbml.LIBSBML_OPERATION_SUCCESS:
                return
            else:
                err_msg = 'Error encountered trying to ' + message + '.' \
                        + 'LibSBML returned error code ' + str(value) + ': "' \
                        + libsbml.OperationReturnValue_toString(value).strip() + '"'
                self.logger.error(err_msg)
                raise SystemExit(err_msg)
        else:
            #self.logger.info(message)
            return None


    def convertToCobra(self):
        try:
            self.cobraModel = cobra.io.read_sbml_model(self.rpsbml.document.toXMLNode().toXMLString(),
                    use_fbc_package=True)
        except cobra.io.sbml.CobraSBMLError as e:
            self.logger.error(e)
            self.logger.error('Cannot convert the libSBML model to Cobra')


    ##
    #
    #   
    def allObj(self, pathId='rp_pathway'):
        fbc_plugin = self.rpsbml.model.getPlugin('fbc')
        self._checklibSBML(fbc_plugin, 'Getting FBC package')
        #print(self.rpsbml.modelName)
        for objId in [i.getId() for i in fbc_plugin.getListOfObjectives()]:
            #print('----> '+str(objId))
            groups = self.rpsbml.model.getPlugin('groups')
            self._checklibSBML(groups, 'Getting groups plugin')
            #ADD ANNOTATIONS IN GROUPS
            #groups_annotations = groups.getAnnotation()
            rp_pathway = groups.getGroup(pathId)
            self._checklibSBML(rp_pathway, 'Getting RP pathway')
            #find all the species of the reactions in the rp_pathway to return the shadow price
            ''' TO BE DETERMINED IF USED 
            mem = []
            for member in rp_pathway.getListOfMembers():
                reac = self.rpsbml.model.getReaction(member.getIdRef())
                for pro in reac.getListOfProducts():
                    mem.append(pro.getSpecies())
                for rea in reac.getListOfReactants():
                    mem.append(rea.getSpecies())
            '''
            #run the FBA
            #fbc_plugin = rpsbml.model.getPlugin('fbc')
            #TODO: switch this to a idependent function and wrap with error handling functions
            self._checklibSBML(fbc_plugin.setActiveObjectiveId(objId), 'Setting active objective '+str(objId))
            self.convertToCobra()
            res = self.cobraModel.optimize()
            #update the annotations with the reaction flux results
            #reac.getChild('RDF').getChild('Ibisba').getChild('ibisba').getChild(child).addAttr('value', str(value))
            #update the annotations with the species shadow prices
            for member in rp_pathway.getListOfMembers():
                reac = self.rpsbml.model.getReaction(member.getIdRef())
                reac_annot = reac.getAnnotation()
                ibisba_annot = reac_annot.getChild('RDF').getChild('Ibisba').getChild('ibisba')
                #NOTE: for some reaseon one neesds to wrap a child XMLnode into a parent node for it to be valid
                #here we just pass the child node to be added
                tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<ibisba:ibisba xmlns:ibisba="http://ibisba.eu"> <ibisba:fba_'+str(objId)+' units="mmol_per_gDW_per_hr" value="'+str(res.fluxes.get(member.getIdRef()))+'" /> </ibisba:ibisba>')
                ibisba_annot.addChild(tmpAnnot.getChild('fba_'+str(objId)))
                ### if targetSink add it to the annotations of the group as well ##
                #ADD ANNOTATIONS IN GROUPS
                #if member.getIdRef()=='targetSink':
                #    groups_annotations.addChild(tmpAnnot.getChild('fba_'+str(objId)))
            del res
            ''' TO BE DETERMINED IF USED
            #update the shadow prices for species
            for speName in list(set(mem)): #remoce duplicates
                #only choose the heterologous species
                if len([x for x in speName.split('_') if x])==4:
                    spe = self.rpsbml.model.getSpecies(speName)
                    spe_annot = spe.getAnnotation()
                    ibisba_annot = spe_annot.getChild('RDF').getChild('Ibisba').getChild('ibisba')
                    tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<ibisba:ibisba xmlns:ibisba="http://ibisba.eu"> <ibisba:shadow_price_'+str(objId)+' units="mmol_per_gDW_per_hr" value="'+str(res.shadow_prices.get(speName))+'" /> </ibisba:ibisba>')
                    ibisba_annot.addChild(tmpAnnot.getChild('shadow_price_'+str(objId)))
            '''


    ## function to add the splitObj or any other objective that does not exist in the model
    #
    #
    def addObjective():
        pass       


    ## Given that there can be multiple objectives defined in the SBML, this function switches between different ones
    #
    # TODO: 
    def switchObjective(self, objId):
        fbc_plugin = self.rpsbml.model.getPlugin('fbc')
        listObj = fbc_plugin.getListOfObjectives()
        if objId in [i.getId() for i in listObj]:
            listObj.setActiveObjective(objId)
            return objId
        else:
            logger.warning('The objective Id '+str(objId)+' does not exist in rpsbml')
            return None
        self.libsbml_to_cobra()


    ##
    #
    #
    def libsbml_to_cobra(self):
        self.cobraModel = cobra.io.read_sbml_model(document.toXMLNode().toXMLString())
        #for an old version of cobrapy (0.4)
        #return cobra.io.sbml3.parse_xml_into_model(lxml.etree.fromstring(rpsbml.document.toXMLNode().toXMLString()))
        
    
    ##
    #
    #TODO: should update the 
    def simulate(self):
        res = self.cobraModel.optimize()
        return self.cobraModel.summary()
        
    
    ##
    #
    # child: rule_score, fba_biomass_score, fba_target_score, fba_splitObj_score, global_score
    def updateFBAscore(self, child, value):
        self.rpsbml()
        #reaction
        for reac in self.rpsbml.getListOfReactions():
            reac.getChild('RDF').getChild('Ibisba').getChild('ibisba').getChild(child).addAttr('value', str(value))
        #pathway


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

#This is to test the good functioing of the main function and make sure that
#there are no errors
#TODO: change this to the stepTest
if __name__ == "__main__":
    #read the TAR.XZ with all the SBML pathways
    rpsbml_paths = {}
    tar = tarfile.open('tests/testFBAin.tar.xz') #TODO: create this
    rpsbml_paths = {}
    for member in tar.getmembers():
        rpsbml_paths[member.name] = rpFBA.rpSBML(member.name,libsbml.readSBMLFromString(tar.extractfile(member).read().decode("utf-8")))
    #pass the different models to the SBML solvers and write the results to file
    #rpsbml_paths = readrpSBMLzip(params.inSBMLzip)
    for rpsbml_name in rpsbml_paths:
        rpfba = rpFBA.rpFBA(rpsbml_paths[rpsbml_name])
        rpfba.allObj()
    #writerpSBMLzip(rpsbml_paths, params.outSBMLzip)
    #designed to write using TAR.XZ with all the SBML pathways
    with tarfile.open('testFBAout.tar.xz', 'w:xz') as tf:
        for rpsbml_name in rpsbml_paths:
            data = libsbml.writeSBMLToString(rpsbml_paths[rpsbml_name].document).encode('utf-8')
            fiOut = BytesIO(data)
            info = tarfile.TarInfo(rpsbml_name)
            info.size = len(data)
            tf.addfile(tarinfo=info, fileobj=fiOut)

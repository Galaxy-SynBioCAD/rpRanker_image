import libsbml
from hashlib import md5
import os
import logging

## @package RetroPath SBML writer
# Documentation for SBML representation of the different model
#
# To exchange between the different workflow nodes, the SBML (XML) format is used. This
# implies using the libSBML library to create the standard definitions of species, reactions, etc...
# Here we also define our own annotations that are used internally in that we call IBISBA nodes.
# The object holds an SBML object and a series of methods to write and access IBISBA related annotations


##################################################################
############################### rpSBML ###########################
##################################################################


## libSBML reader for RetroPath
# Converts an SBML object (or file) into the internal format
#
class rpSBML:
    ## Constructor
    #
    # @param model libSBML model object
    # @param docModel libSBML Document object
    # @param nameSpaceModel libSBML name space (not required)
    def __init__(self, modelName, document=None, path=None):
        self.logger = logging.getLogger(__name__)
        self.logger.info('Started instance of rpSBML')
        self.modelName = modelName
        self.document = document
        if self.document==None:
            self.model = None
        else:
            self.model = self.document.getModel()
        self.path = path
        #need to scan for these if we are passing mode and documents
        self.hetero_group = None
        self.compartmentName = None
        self.compartmentId = None
        #if self.model==None:
        #    self.logger.warning('rpSBML object was initiated as empty. Please call createModel() to initalise the model')

    #######################################################################
    ############################# PRIVATE FUNCTIONS ####################### 
    #######################################################################


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


    ## String to SBML ID
    #
    # Convert any String to one that is compatible with the SBML metaID formatting requirements
    #
    # @param name The input string
    def _nameToSbmlId(self, name):
        IdStream = []
        count = 0
        end = len(name)
        if '0' <= name[count] and name[count] <= '9':
            IdStream.append('_')
        for count in range(0, end):
            if (('0' <= name[count] and name[count] <= '9') or
                    ('a' <= name[count] and name[count] <= 'z') or
                    ('A' <= name[count] and name[count] <= 'Z')):
                IdStream.append(name[count])
            else:
                IdStream.append('_')
        Id = ''.join(IdStream)
        if Id[len(Id) - 1] != '_':
            return Id
        return Id[:-1]


    ## String to hashed ID
    #
    # Hash an input string and then pass it to _nameToSbmlId()
    #
    # @param input string
    def _genMetaID(self, name):
        return self._nameToSbmlId(md5(str(name).encode('utf-8')).hexdigest())


    #####################################################################
    ########################## INPUT/OUTPUT #############################
    #####################################################################


    ## Open an SBML using libSBML 
    #
    # Situation where an SBML is passed to add the heterologous pathway
    #
    # @param inFile String Path to the input SBML file
    def readSBML(self, inFile):
        if not os.path.isfile(inFile):
            self.logger.error('Invalid input file')
            raise FileNotFoundError
        document = libsbml.readSBML(inFile)
        self._checklibSBML(document, 'reading input file')
        errors = document.getNumErrors()
        #display the errors in the log accordning to the severity
        for err in [document.getError(i) for i in range(document.getNumErrors())]:
            if err.isFatal:
                self.logger.error('libSBML reading error: '+str(err.getShortMessage()))
                raise FileNotFoundError
            else:
                self.logger.warning('libSBML reading warning: '+str(err.getShortMessage()))
        model = document.getModel()
        if not model:
            loging.error('Either the file was not read correctly or the SBML is empty')
            raise FileNotFoundError
        self.document = document
        self.model = model


    ## Export a libSBML model to file
    #
    # Export the libSBML model to an SBML file
    # 
    # @param model libSBML model to be saved to file
    # @param model_id model id, note that the name of the file will be that
    # @param path Non required parameter that will define the path where the model will be saved
    def writeSBML(self, path):
        ####### check the path #########
        #need to determine where are the path id's coming from
        p = None
        if path:
            if path[-1:]=='/':
                path = path[:-1]
            if not os.path.isdir(path):
                if self.path:
                    p = self.path
                else:
                    self.logger.error('The output path is not a directory: '+str(path))
                    return False
            else:
                p = path
        else:
            p = self.path
        ########## check and create folder #####
        if not os.path.exists(p):
            os.makedirs(p)
        libsbml.writeSBMLToFile(self.document, p+'/'+str(self.modelName)+'.sbml')
        return True


    #####################################################################
    ########################## READ/WRITE ###############################
    #####################################################################

    
    ## Return the reaction ID's and the pathway annotation
    #
    #
    def readRPpathway(self, pathId='rp_pathway'):
        groups = self.model.getPlugin('groups')
        rp_pathway = groups.getGroup(pathId)
        self._checklibSBML(rp_pathway, 'retreiving groups rp_pathway')
        toRet = []
        for member in rp_pathway.getListOfMembers():
            toRet.append(member.getIdRef())
        return toRet
 

    ## Read the reaction rules from the IBISBA annotation
    #
    #@param path_id default='rp_pathway' unique ID (per SBML) where the heterologous pathways are stored
    #@return toRet dictionnary with the reaction rule and rule_id as key
    def readRPrules(self, path_id='rp_pathway'):
        toRet = {}
        for reacId in self.readRPpathway(path_id):
            reac = self.model.getReaction(reacId)
            ibibsa_annot = self.readIBISBAAnnotation(reac.getAnnotation())
            if not ibibsa_annot['rule_id']=='' and not ibibsa_annot['smiles']=='':
                toRet[ibibsa_annot['rule_id']] = ibibsa_annot['smiles'].replace('&gt;', '>')
        return toRet


    ## Return the species annitations 
    #
    #
    def readRPspecies(self, path_id='rp_pathway'):
        reacMembers = {}
        for reacId in self.readRPpathway(path_id):
            reacMembers[reacId] = {}
            reacMembers[reacId]['products'] = {}
            reacMembers[reacId]['reactants'] = {}
            reac = self.model.getReaction(reacId)
            for pro in reac.getListOfProducts():
                reacMembers[reacId]['products'][pro.getSpecies()] = pro.getStoichiometry()
            for rea in reac.getListOfReactants():
                reacMembers[reacId]['reactants'][rea.getSpecies()] = rea.getStoichiometry()
        return reacMembers


    ## Return the species
    #
    #
    def readUniqueRPspecies(self, pathId='rp_pathway'):
        rpSpecies = self.readRPspecies()
        toRet = []
        for i in rpSpecies:
            for y in rpSpecies[i]:
                for z in rpSpecies[i][y]:
                    if not z in toRet:
                        toRet.append(z)
        return toRet
        #reacMembers = self.readRPspecies(path_id)
        #return set(set(ori_rp_path['products'].keys())|set(ori_rp_path['reactants'].keys()))


    ## Return the MIRIAM annotations of species
    #
    #
    def readMIRIAMAnnotation(self, annot):
        toRet = {}
        bag = annot.getChild('RDF').getChild('Description').getChild('is').getChild('Bag')
        for i in range(bag.getNumChildren()):
            str_annot = bag.getChild(i).getAttrValue(0)
            if str_annot=='':
                self.logger.warning('This contains no attributes: '+str(bag.getChild(i).toXMLString()))
                continue
            dbid = str_annot.split('/')[-2].split('.')[0]
            if len(str_annot.split('/')[-1].split(':'))==2:
                cid = str_annot.split('/')[-1].split(':')[1]
            else:
                cid = str_annot.split('/')[-1]
            if not dbid in toRet: 
                toRet[dbid] = []
            toRet[dbid].append(cid)
        return toRet


    ## Takes for input a libSBML annotatio object and returns a dictionnary of the annotations
    #
    def readIBISBAAnnotation(self, annot):
        toRet = {'dfG_prime_m': {},
                 'dfG_uncert': {},
                 'dfG_prime_o': {},
                 'fba_rpFBA_obj': {},
                 'path_id': None,
                 'step_id': None,
                 'sub_step_id': None,
                 'rule_score': None,
                 'smiles': None,
                 'selenzyme': None,
                 'rule_id': None,
                 'rule_mnxr': None,
                 'rule_score': None,
                 'global_score': None
                }
        bag = annot.getChild('RDF').getChild('Ibisba').getChild('ibisba')
        for i in range(bag.getNumChildren()):
            ann = bag.getChild(i)
            if ann=='':
                self.logger.warning('This contains no attributes: '+str(ann.toXMLString()))
                continue
            #if not ann.getName() in toRet:
            if ann.getName()=='dfG_prime_m' or ann.getName()=='dfG_uncert' or ann.getName()=='dfG_prime_o' or ann.getName()[0:4]=='fba_':
                toRet[ann.getName()] = {
                        'units': ann.getAttrValue('units'), 
                        'value': float(ann.getAttrValue('value'))}
            elif ann.getName()=='path_id' or ann.getName()=='step_id' or ann.getName()=='sub_step_id':
                try:
                    toRet[ann.getName()] = int(ann.getAttrValue('value'))
                except ValueError:
                    toRet[ann.getName()] = None
            elif ann.getName()=='rule_score' or ann.getName()=='global_score':
                try:
                    toRet[ann.getName()] = float(ann.getAttrValue('value'))
                except ValueError:
                    toRet[ann.getName()] = None
            elif ann.getName()=='smiles':
                toRet[ann.getName()] = ann.getChild(0).toXMLString().replace('&gt;', '>')
            elif ann.getName()=='selenzyme':
                toRet['selenzyme'] = {}
                for y in range(ann.getNumChildren()):
                    selAnn = ann.getChild(y)
                    try:
                        toRet['selenzyme'][selAnn.getName()] = float(selAnn.getAttrValue('value'))
                    except ValueError:
                        toRet['selenzyme'][selAnn.getName()] = None
            else:
                toRet[ann.getName()] = ann.getChild(0).toXMLString()
        return toRet


    ## Function to return the products and the species associated with a reaction
    #
    # @return Dictionnary with right==product and left==reactants
    def readReactionSpecies(self, reaction, isID=False):
        #TODO: check that reaction is either an sbml species; if not check that its a string and that
        # it exists in the rpsbml model
        toRet = {'left': {}, 'right': {}}
        #reactants
        for i in range(reaction.getNumReactants()):
            reactant_ref = reaction.getReactant(i)
            reactant = self.model.getSpecies(reactant_ref.getSpecies())
            if isID:
                toRet['left'][reactant.getId()] = int(reactant_ref.getStoichiometry())
            else:
                toRet['left'][reactant.getName()] = int(reactant_ref.getStoichiometry())
        #products
        for i in range(reaction.getNumProducts()):
            product_ref = reaction.getProduct(i)
            product = self.model.getSpecies(product_ref.getSpecies())
            if isID:
                toRet['right'][product.getId()] = int(product_ref.getStoichiometry())
            else:
                toRet['right'][product.getName()] = int(product_ref.getStoichiometry())
            toRet['reversible'] = reaction.getReversible()
        return toRet


    #####################################################################
    ######################### INQUIRE ###################################
    #####################################################################


    ## Function to find out if the model already contains a species according to its name
    #
    #
    def speciesExists(self, speciesName):
        if speciesName in [i.getName() for i in self.model.getListOfSpecies()]:
            return True
        return False


    #########################################################################
    ################### CONVERT BETWEEEN FORMATS ############################
    #########################################################################
    

    ## Really used to complete the monocomponent reactions   
    #{'rule_id': 'RR-01-503dbb54cf91-49-F', 'right': {'TARGET_0000000001': 1}, 'left': {'MNXM2': 1, 'MNXM376': 1}, 'path_id': 1, 'step': 1, 'sub_step': 1, 'transformation_id': 'TRS_0_0_17'}
    #
    def outPathsDict(self, pathId='rp_pathway'):
        pathway = {}
        for member in self.readRPpathway(pathId):
            #TODO: need to find a better way
            if not member=='targetSink':
                reaction = self.model.getReaction(member)
                ibisbaAnnot = self.readIBISBAAnnotation(reaction.getAnnotation())
                speciesReac = self.readReactionSpecies(reaction)
                step = {'reaction_id': member,
                        'reaction_rule': ibisbaAnnot['smiles'], 
                        'rule_score': ibisbaAnnot['rule_score'],
                        'rule_id': ibisbaAnnot['rule_id'],
                        'rule_mnxr': ibisbaAnnot['rule_mnxr'],
                        'right': speciesReac['right'],
                        'left': speciesReac['left'],
                        'path_id': ibisbaAnnot['path_id'],
                        'step': ibisbaAnnot['step_id'],
                        'sub_step': ibisbaAnnot['sub_step_id']}
                pathway[ibisbaAnnot['step_id']] = step
        return pathway


    #########################################################################
    ############################# COMPARE MODELS ############################
    #########################################################################


    ## Find out if two libSBML Species or Reactions come from the same species
    #
    # Compare two dictionnaries and if any of the values of any of the same keys are the same then the 
    # function return True, and if none are found then return False
    #
    # @param libSBML Annotation object for one of the 
    # @return Boolean to determine if they are the same
    def compareIBISBAAnnotations(self, source_annot, target_annot):
        source_dict = self.readIBISBAAnnotation(source_annot)
        target_dict = self.readIBISBAAnnotation(target_annot)
        #list the common keys between the two
        for same_key in list(set(list(source_dict.keys())).intersection(list(target_dict.keys()))):
            if source_dict[same_key]==target_dict[same_key]:
                return True
        return False


    ## Find out if two libSBML Species or Reactions come from the same species
    #
    # Compare two dictionnaries and if any of the values of any of the same keys are the same then the 
    # function return True, and if none are found then return False
    #
    # @param libSBML Annotation object for one of the 
    # @return Boolean to determine if they are the same
    def compareMIRIAMAnnotations(self, source_annot, target_annot):
        source_dict = self.readMIRIAMAnnotation(source_annot)
        target_dict = self.readMIRIAMAnnotation(target_annot)
        #list the common keys between the two
        for com_key in set(list(source_dict.keys()))-(set(list(source_dict.keys()))-set(list(target_dict.keys()))):
            #compare the keys and if same is non-empty means that there 
            #are at least one instance of the key that is the same
            if bool(set(source_dict[com_key]) & set(target_dict[com_key])):
                return True
        return False


    ## Compare an annotation and a dictionnary structured
    #
    #
    def compareAnnotations_annot_dict(self, source_annot, target_dict):
        source_dict = self.readMIRIAMAnnotation(source_annot)
        #list the common keys between the two
        for com_key in set(list(source_dict.keys()))-(set(list(source_dict.keys()))-set(list(target_dict.keys()))):
            #compare the keys and if same is non-empty means that there 
            #are at least one instance of the key that is the same
            if bool(set(source_dict[com_key]) & set(target_dict[com_key])):
                return True
        return False

    
    ## Compare two dictionnaries sutructured as dict
    #
    #
    def compareAnnotations_dict_dict(self, source_dict, target_dict):
        #list the common keys between the two
        for com_key in set(list(source_dict.keys()))-(set(list(source_dict.keys()))-set(list(target_dict.keys()))):
            #compare the keys and if same is non-empty means that there 
            #are at least one instance of the key that is the same
            if bool(set(source_dict[com_key]) & set(target_dict[com_key])):
                return True
        return False


    ## Function to compare two SBML's RP pathways
    #
    # Function that compares the annotations of reactions and if not found, the annotations of all
    # species in that reaction to try to recover the correct ones. Because we are working with
    # intermediate cofactors for the RP generated pathways, the annotation crossreference will
    # not work. Best is to use the cross-reference to the original reaction
    #
    def compareRPpathways(self, measured_sbml):
        #return all the species annotations of the RP pathways
        try:
            meas_rp_species = measured_sbml.readRPspecies()
            found_meas_rp_species = measured_sbml.readRPspecies()
            for meas_step_id in meas_rp_species:
                meas_rp_species[meas_step_id]['annotation'] = measured_sbml.model.getReaction(meas_step_id).getAnnotation()
                found_meas_rp_species[meas_step_id]['found'] = False
                for spe_name in meas_rp_species[meas_step_id]['reactants']:
                    meas_rp_species[meas_step_id]['reactants'][spe_name] = measured_sbml.model.getSpecies(spe_name).getAnnotation()
                    found_meas_rp_species[meas_step_id]['reactants'][spe_name] = False
                for spe_name in meas_rp_species[meas_step_id]['products']:
                    meas_rp_species[meas_step_id]['products'][spe_name] = measured_sbml.model.getSpecies(spe_name).getAnnotation()
                    found_meas_rp_species[meas_step_id]['products'][spe_name] = False
            rp_rp_species = self.readRPspecies()
            for rp_step_id in rp_rp_species:
                rp_rp_species[rp_step_id]['annotation'] = self.model.getReaction(rp_step_id).getAnnotation()
                for spe_name in rp_rp_species[rp_step_id]['reactants']:
                    rp_rp_species[rp_step_id]['reactants'][spe_name] = self.model.getSpecies(spe_name).getAnnotation()
                for spe_name in rp_rp_species[rp_step_id]['products']:
                    rp_rp_species[rp_step_id]['products'][spe_name] = self.model.getSpecies(spe_name).getAnnotation()
            try:
                del rp_rp_species['targetSink']
            except KeyError:
                self.logger.error('The generated RP pathway does not have a targetSink')
        except AttributeError:
            self.logger.error('TODO: debug, for some reason some are passed as None here')
            return False, {}
        #compare the number of steps in the pathway
        if not len(meas_rp_species)==len(rp_rp_species): #add one for the targetSink
            self.logger.warning('The pathways are not of the same length')
            #self.logger.error(len(meas_rp_species))
            #self.logger.error(meas_rp_species.keys())
            #self.logger.error(len(rp_rp_species))
            #self.logger.error(rp_rp_species.keys())
            #self.logger.error(found_meas_rp_species)
            return False, {}
        ############## compare using the reactions ###################
        for meas_step_id in measured_sbml.readRPpathway():
            for rp_step_id in rp_rp_species:
                if self.compareMIRIAMAnnotations(rp_rp_species[rp_step_id]['annotation'], meas_rp_species[meas_step_id]['annotation']):
                    found_meas_rp_species[meas_step_id]['found'] = True
                    found_meas_rp_species[meas_step_id]['rp_step_id'] = rp_step_id
                    #print('FOUND USING REACTION')
                    break
        ############## compare using the species ###################
        for meas_step_id in measured_sbml.readRPpathway():
            #if not found_meas_rp_species[meas_step_id]['found']:
            for rp_step_id in rp_rp_species:
                # We test to see if the meas reaction elements all exist in rp reaction and not the opposite
                #because the measured pathways may not contain all the elements
                ########## reactants ##########
                for meas_spe_id in meas_rp_species[meas_step_id]['reactants']:
                    for rp_spe_id in rp_rp_species[rp_step_id]['reactants']:
                        if self.compareMIRIAMAnnotations(meas_rp_species[meas_step_id]['reactants'][meas_spe_id], rp_rp_species[rp_step_id]['reactants'][rp_spe_id]):
                            found_meas_rp_species[meas_step_id]['reactants'][meas_spe_id] = True
                            break
                        else:
                            if self.compareIBISBAAnnotations(meas_rp_species[meas_step_id]['reactants'][meas_spe_id], rp_rp_species[rp_step_id]['reactants'][rp_spe_id]):
                                found_meas_rp_species[meas_step_id]['reactants'][meas_spe_id] = True
                                break
                ########### products ###########
                for meas_spe_id in meas_rp_species[meas_step_id]['products']:
                    for rp_spe_id in rp_rp_species[rp_step_id]['products']:
                        if self.compareMIRIAMAnnotations(meas_rp_species[meas_step_id]['products'][meas_spe_id], rp_rp_species[rp_step_id]['products'][rp_spe_id]):
                            found_meas_rp_species[meas_step_id]['products'][meas_spe_id] = True
                            break
                        else:
                            if self.compareIBISBAAnnotations(meas_rp_species[meas_step_id]['products'][meas_spe_id], rp_rp_species[rp_step_id]['products'][rp_spe_id]):
                                found_meas_rp_species[meas_step_id]['products'][meas_spe_id] = True
                                break
                ######### test to see the difference
                pro_found = [found_meas_rp_species[meas_step_id]['products'][i] for i in found_meas_rp_species[meas_step_id]['products']]
                rea_found = [found_meas_rp_species[meas_step_id]['reactants'][i] for i in found_meas_rp_species[meas_step_id]['reactants']]
                if pro_found and rea_found:
                    if all(pro_found) and all(rea_found):
                        found_meas_rp_species[meas_step_id]['found'] = True
                        found_meas_rp_species[meas_step_id]['rp_step_id'] = rp_step_id
                        break
        ################# Now see if all steps have been found ############
        if all(found_meas_rp_species[i]['found'] for i in found_meas_rp_species):
            found_meas_rp_species['measured_model_id'] = measured_sbml.model.getId()
            found_meas_rp_species['rp_model_id'] = self.model.getId()
            return True, found_meas_rp_species
        else:
            return False, {}


    #########################################################################
    ############################# MODEL APPEND ##############################
    #########################################################################


    ## Merge two models species and reactions using the annotations to recognise the same species and reactions
    #
    # The source mode has to have both the GROUPS and FBC packages enabled in its SBML. The course must have a groups
    #called rp_pathway.
    # We add the reactions and species from the rpsbml to the target_model
    #
    # @param bilevel_obj: Defines if the different objectves are to be added to the model (False if (0.0, 0.0) and what are the coefficients to add it)
    #
    def mergeModels(self, target_model, pathId='rp_pathway', bilevel_obj=(0.0, 0.0)):
        #target_model = target_document.getModel()
        #Find the ID's of the similar target_model species
        ################ UNITDEFINITIONS ######
        #return the list of unit definitions id's for the target to avoid overwritting
        #WARNING: this means that the original unit definitions will be prefered over the new one
        target_unitDefID = [i.getId() for i in target_model.getListOfUnitDefinitions()]
        for source_unitDef in self.model.getListOfUnitDefinitions():
            if not source_unitDef.getId() in target_unitDefID: #have to compare by ID since no annotation
                #create a new unitDef in the target
                target_unitDef = target_model.createUnitDefinition()
                self._checklibSBML(target_unitDef, 'fetching target unit definition')
                #copy unitDef info to the target
                self._checklibSBML(target_unitDef.setId(source_unitDef.getId()), 
                    'setting target unit definition ID')
                #self._checklibSBML(target_unitDef.setMetaId(source_unitDef.getMetaId()), 
                #    'setting target unit definition MetaId')
                self._checklibSBML(target_unitDef.setAnnotation(source_unitDef.getAnnotation()), 
                    'setting target unit definition Annotation')
                for source_unit in source_unitDef.getListOfUnits():
                    #copy unit info to the target unitDef
                    target_unit = target_unitDef.createUnit()
                    self._checklibSBML(target_unit, 'creating target unit')
                    self._checklibSBML(target_unit.setKind(source_unit.getKind()), 
                        'setting target unit kind')
                    self._checklibSBML(target_unit.setExponent(source_unit.getExponent()), 
                        'setting target unit exponent')
                    self._checklibSBML(target_unit.setScale(source_unit.getScale()), 
                        'setting target unit scale')
                    self._checklibSBML(target_unit.setMultiplier(source_unit.getMultiplier()), 
                        'setting target unit multiplier')
                target_unitDefID.append(source_unitDef.getId()) #add to the list to make sure its not added twice
        ################ COMPARTMENTS ###############
        sourceCompartmentID_targetCompartmentID = {}
        toAddNum = []
        ####### compare the annotations to find the same ones #######
        for i in range(self.model.getNumCompartments()):
            found = False
            source_compartment = self.model.getCompartment(i)
            self._checklibSBML(source_compartment, 'Getting target compartment')
            source_annotation = source_compartment.getAnnotation()
            #self._checklibSBML(source_annotation, 'Getting compartment target annotation')
            if not source_annotation:
                self.logger.warning('No annotation for the source of compartment '+str(source_compartment.getId()))
                continue
            for y in range(target_model.getNumCompartments()):
                target_compartment = target_model.getCompartment(y)
                self._checklibSBML(target_compartment, 'Getting target compartment')
                target_annotation = target_compartment.getAnnotation()
                #self._checklibSBML(target_annotation, 'Getting target annotation')
                if not target_annotation:
                    self.logger.warning('No annotation for the target of compartment: '+str(target_compartment.getId()))
                    continue
                if self.compareMIRIAMAnnotations(source_annotation, target_annotation):
                    sourceCompartmentID_targetCompartmentID[source_compartment.getId()] = target_compartment.getId() 
                    found = True
                    break
            if not found:
                toAddNum.append(i)
        for i in toAddNum:
            source_compartment = self.model.getCompartment(i)
            self._checklibSBML(source_compartment, 'Getting target compartment')
            target_compartment = target_model.createCompartment()
            self._checklibSBML(target_compartment, 'Creating target compartment')
            self._checklibSBML(target_compartment.setMetaId(source_compartment.getMetaId()),
                    'setting target metaId')
            self._checklibSBML(target_compartment.setId(source_compartment.getId()),
                    'setting target id')
            self._checklibSBML(target_compartment.setName(source_compartment.getName()),
                    'setting target name')
            self._checklibSBML(target_compartment.setConstant(source_compartment.getConstant()),
                    'setting target constant')
            self._checklibSBML(target_compartment.setAnnotation(source_compartment.getAnnotation()),
                    'setting target annotation')
            self._checklibSBML(target_compartment.setSBOTerm(source_compartment.getSBOTerm()),
                    'setting target annotation')
        ################ PARAMETERS ###########
        #WARNING: here we compare by ID
        targetParametersID = [i.getId() for i in target_model.getListOfParameters()]
        for source_parameter in self.model.getListOfParameters():
            if not source_parameter.getId() in targetParametersID:
                target_parameter = target_model.createParameter()
                self._checklibSBML(target_parameter, 'creating target parameter')
                self._checklibSBML(target_parameter.setId(source_parameter.getId()), 'setting target parameter ID')
                self._checklibSBML(target_parameter.setSBOTerm(source_parameter.getSBOTerm()), 
                    'setting target parameter SBO')
                self._checklibSBML(target_parameter.setUnits(source_parameter.getUnits()),
                    'setting target parameter Units')
                self._checklibSBML(target_parameter.setValue(source_parameter.getValue()),
                    'setting target parameter Value') 
                self._checklibSBML(target_parameter.setConstant(source_parameter.getConstant()),
                    'setting target parameter ID')
        ################ MODEL FBC ########################
        if not target_model.isPackageEnabled('fbc'):
            self._checklibSBML(target_model.enablePackage(
                'http://www.sbml.org/sbml/level3/version1/fbc/version2',
                'fbc',
                True),
                    'Enabling the FBC package')
        target_fbc = target_model.getPlugin('fbc')
        source_fbc = self.model.getPlugin('fbc')
        ################ FBC GENE PRODUCTS ########################
        #make a list of all the gene product
        #WARNING: here we compare by ID
        targetGenProductID = [i.getId() for i in target_fbc.getListOfGeneProducts()]
        for source_geneProduct in source_fbc.getListOfGeneProducts():
            if not source_geneProduct.getId() in targetGenProductID:
                target_geneProduct = target_fbc.createGeneProduct()
                self._checklibSBML(target_geneProduct, 'creating target gene product')
                self._checklibSBML(target_geneProduct.setId(source_geneProduct.getId()), 
                    'setting target gene product id')
                self._checklibSBML(target_geneProduct.setLabel(source_geneProduct.getLabel()), 
                    'setting target gene product label')
                self._checklibSBML(target_geneProduct.setName(source_geneProduct.getName()), 
                    'setting target gene product name')
                self._checklibSBML(target_geneProduct.setMetaId(source_geneProduct.getMetaId()),
                    'setting target gene product metaID')
        ############### FBC OBJECTIVES ############
        #WARNING: here we compare the Objective by ID, and we add the downstream fluxObjectives
        targetObjectiveID = [i.getId() for i in target_fbc.getListOfObjectives()]
        sourceObjectiveID = [i.getId() for i in source_fbc.getListOfObjectives()]
        for source_objective in source_fbc.getListOfObjectives():
            if not source_objective.getId() in targetObjectiveID:
                target_objective = target_fbc.createObjective()
                self._checklibSBML(target_objective, 'creating target objective')
                self._checklibSBML(target_objective.setId(source_objective.getId()), 'setting target objective')
                self._checklibSBML(target_objective.setName(source_objective.getName()), 'setting target objective')
                self._checklibSBML(target_objective.setType(source_objective.getType()), 
                        'setting target objective type')
                for source_fluxObjective in source_objective.getListOfFluxObjectives():
                    target_fluxObjective = target_objective.createFluxObjective()
                    self._checklibSBML(target_fluxObjective, 'creating target flux objective')
                    self._checklibSBML(target_fluxObjective.setName(source_fluxObjective.getName()),
                        'setting target flux objective name')
                    self._checklibSBML(target_fluxObjective.setCoefficient(source_fluxObjective.getCoefficient()),
                        'setting target flux objective coefficient')
                    self._checklibSBML(target_fluxObjective.setReaction(source_fluxObjective.getReaction()),
                        'setting target flux objective reaction')
        #### bilevel
        #add a bilevel fluxObjective. Under the assumption that the GEM model input only has a single objective and that that is the biomass one. 
        #NOTE: only if each model have only one objective with a single flux objective
        #create the two list of objectives
        if 1>=len(targetObjectiveID)>0 and 1>=len(sourceObjectiveID)>0:
            #create a new one objective
            bilevel_objective = target_fbc.createObjective()
            self._checklibSBML(bilevel_objective, 'creating bilevel objective')
            self._checklibSBML(bilevel_objective.setId('rpFBA_bilevel_obj'), 'setting bilevel obj id')
            self._checklibSBML(bilevel_objective.setType('maximize'), 'setting type of bilevel')
            #list the flux objectives
            source_fluxObjective = source_fbc.getObjective(sourceObjectiveID[0])
            target_fluxObjective = target_fbc.getObjective(targetObjectiveID[0])
            targetFluxObjectives = target_fluxObjective.getListOfFluxObjectives()
            sourceFluxObjectives = source_fluxObjective.getListOfFluxObjectives()
            if 1>=len(targetFluxObjectives)>0 and 1>=len(sourceFluxObjectives)>0:
                #biomass
                bilevel_fluxObjective_biomass = bilevel_objective.createFluxObjective()
                self._checklibSBML(bilevel_fluxObjective_biomass.setName(sourceFluxObjectives[0].getName()),
                    'setting target flux objective name')
                self._checklibSBML(bilevel_fluxObjective_biomass.setCoefficient(0.5),
                    'setting target flux objective coefficient')
                self._checklibSBML(bilevel_fluxObjective_biomass.setReaction(sourceFluxObjectives[0].getReaction()),
                    'setting target flux objective reaction')
                #target
                bilevel_fluxObjective_target = bilevel_objective.createFluxObjective()
                self._checklibSBML(bilevel_fluxObjective_target, 'creating target flux objective')
                self._checklibSBML(bilevel_fluxObjective_target.setName(targetFluxObjectives[0].getName()),
                    'setting target flux objective name')
                self._checklibSBML(bilevel_fluxObjective_target.setCoefficient(0.5),
                    'setting target flux objective coefficient')
                self._checklibSBML(bilevel_fluxObjective_target.setReaction(targetFluxObjectives[0].getReaction()),
                    'setting target flux objective reaction')
            else:
                self.logger.warning('Either the target or source model has one of the oobjectives with multiple flux values')
        else:
            self.logger.warning('There are more than one, or zero objective in the target and the source')
        ################ SPECIES ####################
        #TODO: modify the name to add rpPathway
        sourceSpeciesID_targetSpeciesID = {}
        #### compare
        # first make the target model dictionnary of the species for the target model
        targetModel_speciesAnnot = {}
        for y in range(target_model.getNumSpecies()):
            #target_species = target_model.getSpecies(y)
            #self._checklibSBML(target_species, 'Getting target species')
            #target_annotation = target_species.getAnnotation()
            target_annotation = target_model.getSpecies(y).getAnnotation()
            if not target_annotation:    
                self.logger.warning('Cannot find annotations for species: '+str(target_model.getSpecies(y).getId()))
                continue
            self._checklibSBML(target_annotation, 'Getting target annotation')
            targetModel_speciesAnnot[y] = self.readMIRIAMAnnotation(target_annotation)
        # second make the target model dictionnary for the species of the rp model
        toAddNum = []
        for i in range(self.model.getNumSpecies()):
            found = False
            source_species = self.model.getSpecies(i)
            self._checklibSBML(source_species, 'Getting source species')
            source_annotation = source_species.getAnnotation()
            self._checklibSBML(source_annotation, 'Getting source annotation')
            if not source_annotation:
                self.logger.warning('No annotation for the source of compartment '+str(source_compartment.getId()))
                #we assume that if there are no annotations then we add it
                toAddNum.append(i)
                continue
            #for y in range(target_model.getNumSpecies()):
                #target_species = target_model.getSpecies(y)
                #self._checklibSBML(target_species, 'Getting target species')
                #target_annotation = target_species.getAnnotation()
                #self._checklibSBML(target_annotation, 'Getting target annotation')
                #if not target_annotation:    
                #    self.logger.warning('Cannot find target number: '+str(y))
                #    continue
            for y in targetModel_speciesAnnot:
                #if self.compareMIRIAMAnnotations(source_annotation, targetModel_speciesAnnot[y]):
                if self.compareAnnotations_annot_dict(source_annotation, targetModel_speciesAnnot[y]):
                    #save the speciesID as being the same
                    #sourceSpeciesID_targetSpeciesID[self.model.species[i].getId()] = target_model.species[y].getId()
                    sourceSpeciesID_targetSpeciesID[self.model.species[i].getId()] = target_model.getSpecies(y).getId()
                    found = True
                    break
            #if it has not been found then add it to the target_model
            if not found:
                toAddNum.append(i)
        for i in toAddNum:
            source_species = self.model.getSpecies(i)
            self._checklibSBML(source_species, 'fetching source species')
            target_species = target_model.createSpecies()
            self._checklibSBML(target_species, 'creating species')
            self._checklibSBML(target_species.setMetaId(source_species.getMetaId()),
                    'setting target metaId')
            self._checklibSBML(target_species.setId(source_species.getId()),
                    'setting target id')
            self._checklibSBML(target_species.setCompartment(source_species.getCompartment()),
                    'setting target compartment')
            self._checklibSBML(target_species.setInitialConcentration(
                source_species.getInitialConcentration()),
                    'setting target initial concentration')
            self._checklibSBML(target_species.setBoundaryCondition(
                source_species.getBoundaryCondition()),
                    'setting target boundary concentration')
            self._checklibSBML(target_species.setHasOnlySubstanceUnits(
                source_species.getHasOnlySubstanceUnits()),
                    'setting target has only substance units')
            self._checklibSBML(target_species.setBoundaryCondition(
                source_species.getBoundaryCondition()),
                    'setting target boundary condition')
            self._checklibSBML(target_species.setConstant(source_species.getConstant()),
                'setting target constant')
            self._checklibSBML(target_species.setAnnotation(source_species.getAnnotation()),
                'setting target annotation')
        ################ REACTIONS ###################
        #Find the ID's of the similar target_model reactions
        #need to create a new instance of reactions to add to the model
        #test to see if the target model has the FBC package and if not add it
        sourceReactionsID_targetReactionsID = {}
        toAddNum = []
        if not target_model.isPackageEnabled('fbc'):
            self._checklibSBML(target_model.enablePackage(
                'http://www.sbml.org/sbml/level3/version1/fbc/version2',
                'fbc',
                True),
                    'Enabling the FBC package')
        #note sure why one needs to set this as False
        self._checklibSBML(self.document.setPackageRequired('fbc', False), 'enabling FBC package')
        #### compare the annotations to find the co-factors
        #first make the target model dictionnary of the reactions
        targetModel_reactionsAnnot = {}
        for y in range(target_model.getNumReactions()):
            #target_reaction = target_model.getReaction(y)
            #self._checklibSBML(target_reaction, 'fetching target reaction annotation') 
            #target_annotation = target_reaction.getAnnotation()
            #self._checklibSBML(target_annotation, 'fetching target reaction annotation') 
            target_annotation = target_model.getReaction(y).getAnnotation()
            if not target_annotation:    
                self.logger.warning('No annotation for the target of reaction: '+str(target_model.getReaction(y).getId()))
                continue
            self._checklibSBML(target_annotation, 'fetching target reaction annotation') 
            targetModel_reactionsAnnot[y] = self.readMIRIAMAnnotation(target_annotation)
        #### WANRING: important to list the heterologous pathways in the original model and if
        # comparing the annotations returns true to not add them
        # this is a fix to a bug caused by adding EC numbers to the reactions
        model_rpPathway = self.readRPpathway()
        for i in range(self.model.getNumReactions()):
            found = False
            source_reaction = self.model.getReaction(i)
            self._checklibSBML(source_reaction, 'fetching source reaction')
            source_annotation = source_reaction.getAnnotation()
            #self._checklibSBML(source_annotation, 'fetching source reaction annotation')
            if not source_annotation:
                self.logger.warning(source_annotation)
                self.logger.warning('No annotation for the source of reaction: '+str(source_reaction.getId()))
                toAddNum.append(i)
                continue
            if source_reaction.getId() in model_rpPathway:
                toAddNum.append(i)
                continue
            #for y in range(target_model.getNumReactions()):
                #target_reaction = target_model.getReaction(y)
                #self._checklibSBML(target_reaction, 'fetching target reaction annotation') 
                #target_annotation = target_reaction.getAnnotation()
                #self._checklibSBML(target_annotation, 'fetching target reaction annotation') 
                #if not target_annotation:    
                #    self.logger.warning('No annotation for the target of reaction: '+str(target_reaction.getId()))
                #    continue
            for y in targetModel_reactionsAnnot:
                #if self.compareMIRIAMAnnotations(source_annotation, target_annotation):
                #if self.compareMIRIAMAnnotations(source_annotation, targetModel_reactionsAnnot[y]):
                if self.compareAnnotations_annot_dict(source_annotation, targetModel_reactionsAnnot[y]):
                    #sourceReactionsID_targetReactionsID[self.model.reactions[i].getId()] = target_model.reactions[y].getId()
                    sourceReactionsID_targetReactionsID[source_reaction.getId()] = target_model.getReaction(y).getId()
                    found = True
                    break
            if not found:
                toAddNum.append(i)
        for i in toAddNum:
            source_reaction = self.model.getReaction(i)
            self._checklibSBML(source_reaction, 'fetching source reaction')
            target_reaction = target_model.createReaction()
            self._checklibSBML(target_reaction, 'create reaction')
            target_fbc = target_reaction.getPlugin('fbc')
            self._checklibSBML(target_fbc, 'fetching target FBC package')
            source_fbc = source_reaction.getPlugin('fbc')
            self._checklibSBML(source_fbc, 'fetching source FBC package')
            source_upperFluxBound = source_fbc.getUpperFluxBound()
            self._checklibSBML(source_upperFluxBound, 'fetching upper flux bound')
            self._checklibSBML(target_fbc.setUpperFluxBound(source_upperFluxBound), 
                    'setting upper flux bound')
            source_lowerFluxBound = source_fbc.getLowerFluxBound()
            self._checklibSBML(source_lowerFluxBound, 'fetching lower flux bound')
            self._checklibSBML(target_fbc.setLowerFluxBound(source_lowerFluxBound), 
                    'setting lower flux bound')
            self._checklibSBML(target_reaction.setId(source_reaction.getId()), 'set reaction id')
            self._checklibSBML(target_reaction.setName(source_reaction.getName()), 'set name')
            self._checklibSBML(target_reaction.setSBOTerm(source_reaction.getSBOTerm()), 
                    'setting the reaction system biology ontology (SBO)') #set as process
            #TODO: consider having the two parameters as input to the function
            self._checklibSBML(target_reaction.setReversible(source_reaction.getReversible()), 
                    'set reaction reversibility flag')
            self._checklibSBML(target_reaction.setFast(source_reaction.getFast()), 
                    'set reaction "fast" attribute')
            self._checklibSBML(target_reaction.setMetaId(source_reaction.getMetaId()), 'setting species metaID')
            self._checklibSBML(target_reaction.setAnnotation(source_reaction.getAnnotation()), 
                    'setting annotation for source reaction')
            #reactants_dict
            for y in range(source_reaction.getNumReactants()):
                target_reactant = target_reaction.createReactant()
                self._checklibSBML(target_reactant, 'create target reactant')
                source_reactant = source_reaction.getReactant(y)
                self._checklibSBML(source_reactant, 'fetch source reactant')
                try:
                    #try to get the reactant from the dictionnary if annotations comparison 
                    #elects them to be the same
                    reactantID = sourceSpeciesID_targetSpeciesID[source_reactant.species]
                except KeyError:
                    #if not found in dictionnary then muct be part of the added ones
                    reactantID = source_reactant.species
                self._checklibSBML(target_reactant.setSpecies(reactantID), 'assign reactant species')
                #TODO: check to see the consequences of heterologous parameters not being constant
                self._checklibSBML(target_reactant.setConstant(source_reactant.getConstant()), 
                        'set "constant" on species '+str(source_reactant.getConstant()))
                self._checklibSBML(target_reactant.setStoichiometry(source_reactant.getStoichiometry()),
                        'set stoichiometry ('+str(source_reactant.getStoichiometry)+')')
            #products_dict
            for y in range(source_reaction.getNumProducts()):
                target_product = target_reaction.createProduct()
                self._checklibSBML(target_product, 'create target product')
                source_product = source_reaction.getProduct(y)
                self._checklibSBML(source_product, 'fetch source product')
                try:
                    #try to get the reactant from the dictionnary if annotations comparison 
                    #elects them to be the same
                    productID = sourceSpeciesID_targetSpeciesID[source_product.species]
                except KeyError:
                    #if not found in dictionnary then muct be part of the added ones
                    productID = source_product.species
                self._checklibSBML(target_product.setSpecies(productID), 'assign reactant species')
                #TODO: check to see the consequences of heterologous parameters not being constant
                self._checklibSBML(target_product.setConstant(source_product.getConstant()), 
                        'set "constant" on species '+str(source_product.getConstant()))
                self._checklibSBML(target_product.setStoichiometry(source_product.getStoichiometry()),
                        'set stoichiometry ('+str(source_product.getStoichiometry)+')')
        #### GROUPS #####
        if not target_model.isPackageEnabled('groups'):
            self._checklibSBML(target_model.enablePackage(
                'http://www.sbml.org/sbml/level3/version1/groups/version1',
                'groups',
                True),
                    'Enabling the GROUPS package')
        #!!!! must be set to false for no apparent reason
        self._checklibSBML(self.document.setPackageRequired('groups', False), 'enabling groups package')
        #self._checklibSBML(self.sbmlns.addPkgNamespace('groups',1), 'Add groups package')
        source_groups = self.model.getPlugin('groups')
        self._checklibSBML(source_groups, 'fetching the source model groups')
        target_groups = target_model.getPlugin('groups')
        self._checklibSBML(target_groups, 'fetching the target model groups')
        self._checklibSBML(target_groups.addGroup(source_groups.getGroup(pathId)),
                'copying the source groups "rp_pathway" to the target groups')
        #return the fluxObj for the original model to define the bilevel objective        
        ###### TITLES #####
        target_model.setId(target_model.getId()+'_'+self.model.getId())
        target_model.setName(target_model.getName()+' merged with '+self.model.getId())


    #########################################################################
    ############################# MODEL CREATION FUNCTIONS ##################
    #########################################################################


    ## Create libSBML model instance
    #
    # Function that creates a new libSBML model instance and initiates it with the appropriate packages. Creates a cytosol compartment
    #
    # @param name The name of the model
    # @param modelID The id of the mode
    # @param metaID metaID of the model. Default None means that we will generate a hash from the modelID
    def createModel(self, name, modelID, metaID=None):
        ## sbmldoc
        self.sbmlns = libsbml.SBMLNamespaces(3,1)
        self._checklibSBML(self.sbmlns, 'generating model namespace')
        self._checklibSBML(self.sbmlns.addPkgNamespace('groups',1), 'Add groups package')
        self._checklibSBML(self.sbmlns.addPkgNamespace('fbc',2), 'Add FBC package')
        #sbmlns = libsbml.SBMLNamespaces(3,1,'groups',1)
        self.document = libsbml.SBMLDocument(self.sbmlns)
        self._checklibSBML(self.document, 'generating model doc')
        #!!!! must be set to false for no apparent reason
        self._checklibSBML(self.document.setPackageRequired('fbc', False), 'enabling FBC package') 
        #!!!! must be set to false for no apparent reason
        self._checklibSBML(self.document.setPackageRequired('groups', False), 'enabling groups package')
        ## sbml model
        self.model = self.document.createModel()
        self._checklibSBML(self.model, 'generating the model')
        self._checklibSBML(self.model.setId(modelID), 'setting the model ID')
        model_fbc = self.model.getPlugin('fbc')
        model_fbc.setStrict(True)
        if metaID==None:
            metaID = self._genMetaID(modelID)
        self._checklibSBML(self.model.setMetaId(metaID), 'setting model metaID')
        self._checklibSBML(self.model.setName(name), 'setting model name')
        self._checklibSBML(self.model.setTimeUnits('second'), 'setting model time unit')
        self._checklibSBML(self.model.setExtentUnits('mole'), 'setting model compartment unit')
        self._checklibSBML(self.model.setSubstanceUnits('mole'), 'setting model substance unit')


    ## Create libSBML compartment 
    #
    # cytoplasm compartment TODO: consider seperating it in another function if another compartment is to be created
    #
    # @param model libSBML model object to add the compartment
    # @param size Set the compartement size
    # @return boolean Execution success
    #TODO: set the compName as None by default. To do that you need to regenerate the compXref to 
    #use MNX ids as keys instead of the string names
    def createCompartment(self, size, compId, compName, compXref, metaID=None):
        comp = self.model.createCompartment()
        self._checklibSBML(comp, 'create compartment')
        self._checklibSBML(comp.setId(compId), 'set compartment id')
        self.compartmentId = compId
        if compName:
            self._checklibSBML(comp.setName(compName), 'set the name for the cytoplam')
            self.compartmentName = compName
        self._checklibSBML(comp.setConstant(True), 'set compartment "constant"')
        self._checklibSBML(comp.setSize(size), 'set compartment "size"')
        self._checklibSBML(comp.setSBOTerm(290), 'set SBO term for the cytoplasm compartment')
        if metaID==None:
            metaID = self._genMetaID(compId)
        self._checklibSBML(comp.setMetaId(metaID), 'set the metaID for the compartment')
        #self.compartmentNames.append(compartmentName) #this assumes there is only one compartmentName
        annotation = '''<annotation>
  <rdf:RDF 
  xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" 
  xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" 
  xmlns:bqmodel="http://biomodels.net/model-qualifiers/">'''
        # if the name of the species is MNX then we annotate it using MIRIAM compliance
        #TODO: need to add all known xref from different databases (not just MetaNetX)
        annotation += '''
    <rdf:Description rdf:about="#'''+str(metaID or '')+'''">
      <bqbiol:is>
        <rdf:Bag>'''
        #TODO: for yout to complete
        id_ident = {'mnx': 'metanetx.compartment/', 'bigg': 'bigg.compartment/', 'seed': 'seed/', 'name': 'name/'}
        #WARNING: compartmentNameID as of now, needs to be a MNX ID
        for databaseId in compXref:
            for compartmentId in compXref[databaseId]:
                try:
                    annotation += '''
      <rdf:li rdf:resource="http://identifiers.org/'''+str(id_ident[databaseId])+str(compartmentId)+'''"/>'''
                except KeyError:
                    continue
        annotation += '''
        </rdf:Bag>
      </bqbiol:is>
    </rdf:Description> 
  </rdf:RDF>
</annotation>'''
        self._checklibSBML(comp.setAnnotation(annotation), 'setting annotation for reaction '+str(compName))


    ## Create libSBML unit definition
    #
    # Function that creates a unit definition (composed of one or more units)
    #
    # @param model libSBML model to add the unit definition
    # @param unit_id ID for the unit definition
    # @param metaID metaID for the unit definition. If None creates a hash from unit_id
    # @return Unit definition
    def createUnitDefinition(self, unit_id, metaID=None):
        unitDef = self.model.createUnitDefinition()
        self._checklibSBML(unitDef, 'creating unit definition')
        self._checklibSBML(unitDef.setId(unit_id), 'setting id')
        if metaID==None:
            metaID = self._genMetaID(unit_id)
        self._checklibSBML(unitDef.setMetaId(metaID), 'setting metaID')
        #self.unitDefinitions.append(unit_id)
        return unitDef


    ## Create libSBML unit
    #
    # Function that created a unit
    #
    # @param unitDef libSBML unit definition
    # @param libsmlunit libSBML unit parameter
    # @param exponent Value for the exponent (ex 10^5 mol/sec)
    # @param scale Value for the scale 
    # @param multiplier Value for the multiplie
    # @return Unit
    def createUnit(self, unitDef, libsbmlunit, exponent, scale, multiplier):
        unit = unitDef.createUnit()
        self._checklibSBML(unit, 'creating unit')
        self._checklibSBML(unit.setKind(libsbmlunit), 'setting the kind of unit')
        self._checklibSBML(unit.setExponent(exponent), 'setting the exponenent of the unit')
        self._checklibSBML(unit.setScale(scale), 'setting the scale of the unit')
        self._checklibSBML(unit.setMultiplier(multiplier), 'setting the multiplier of the unit')


    ## Create libSBML parameters
    #
    # Parameters, in our case, are used for the bounds for FBA analysis. Unit parameter must be an instance of unitDefinition
    #
    # @param parameter_id SBML id
    # @param value Float value for this parameter
    # @param unit libSBML unit parameter
    # @param metaID String Optional parameter for SBML metaID
    # @return libSBML parameter object
    def createParameter(self, parameter_id, value, unit, metaID=None):
        newParam = self.model.createParameter()
        self._checklibSBML(newParam, 'Creating a new parameter object')
        self._checklibSBML(newParam.setConstant(True), 'setting as constant')
        self._checklibSBML(newParam.setId(parameter_id), 'setting ID')
        self._checklibSBML(newParam.setValue(value), 'setting value')
        self._checklibSBML(newParam.setUnits(unit), 'setting units')
        self._checklibSBML(newParam.setSBOTerm(625), 'setting SBO term')
        if metaID==None:
            metaID = self._genMetaID(parameter_id)
        self._checklibSBML(newParam.setMetaId(metaID), 'setting meta ID')
        #self.parameters.append(parameter_id)
        return newParam


    ## Create libSBML reaction
    #
    # Create a reaction. fluxBounds is a list of libSBML.UnitDefinition, length of exactly 2 with the first position that is the upper bound and the second is the lower bound. reactants_dict and reactants_dict are dictionnaries that hold the following parameters: name, compartment, stoichiometry
    #
    # @param name Name for the reaction
    # @param reaction_id Reaction ID
    # @param fluxUpperBounds FBC id for the upper flux bound for this reaction
    # @param fluxLowerBounds FBC id for the lower flux bound for this reaction
    # BILAL check the lower
    # @param step 2D dictionnary with the following structure {'left': {'name': stoichiometry, ...}, 'right': {}}
    # @param reaction_smiles String smiles description of this reaction (added in IBISBA annotation)
    # @param compartmentId String Optinal parameter compartment ID
    # @param isTarget Boolean Flag to suppress the warning that the passed step is missing information. Used in this case for the target compound
    # @param hetero_group Groups Optional parameter object that holds all the heterologous pathways
    # @param metaID String Optional parameter reaction metaID
    # @return metaID meta ID for this reaction
    def createReaction(self,
            reacId,
            fluxUpperBound,
            fluxLowerBound,
            step,
            compartmentId,
            reaction_smiles=None,
            ecs=[],
            reacXref={},
            hetero_group=None,
            metaID=None):
        reac = self.model.createReaction()
        self._checklibSBML(reac, 'create reaction')
        ################ FBC ####################
        reac_fbc = reac.getPlugin('fbc')
        self._checklibSBML(reac_fbc, 'extending reaction for FBC')
        #bounds
        self._checklibSBML(reac_fbc.setUpperFluxBound(fluxUpperBound), 'setting '+str(reacId)+' upper flux bound')
        self._checklibSBML(reac_fbc.setLowerFluxBound(fluxLowerBound), 'setting '+str(reacId)+' lower flux bound')
        #########################################
        #reactions
        self._checklibSBML(reac.setId(reacId), 'set reaction id') #same convention as cobrapy
        #self._checklibSBML(reac.setName(str(reacId)+), 'set name') #same convention as cobrapy
        self._checklibSBML(reac.setSBOTerm(185), 'setting the system biology ontology (SBO)') #set as process
        #TODO: consider having the two parameters as input to the function
        self._checklibSBML(reac.setReversible(True), 'set reaction reversibility flag')
        self._checklibSBML(reac.setFast(False), 'set reaction "fast" attribute')
        if metaID==None:
            metaID = self._genMetaID(reacId)
        self._checklibSBML(reac.setMetaId(metaID), 'setting species metaID')
        #reactants_dict
        for reactant in step['left']:
            spe = reac.createReactant()
            self._checklibSBML(spe, 'create reactant')
            #use the same writing convention as CobraPy
            if compartmentId:
                self._checklibSBML(spe.setSpecies(str(reactant)+'__64__'+str(compartmentId)), 'assign reactant species')
            else:
                self._checklibSBML(spe.setSpecies(str(reactant)+'__64__'+str(self.compartmentId)), 'assign reactant species')
            #TODO: check to see the consequences of heterologous parameters not being constant
            self._checklibSBML(spe.setConstant(True), 'set "constant" on species '+str(reactant))
            self._checklibSBML(spe.setStoichiometry(float(step['left'][reactant])),
                'set stoichiometry ('+str(float(step['left'][reactant]))+')')
        #products_dict
        for product in step['right']:
            pro = reac.createProduct()
            self._checklibSBML(pro, 'create product')
            if compartmentId:
                self._checklibSBML(pro.setSpecies(str(product)+'__64__'+str(compartmentId)), 'assign product species')
            else:
                self._checklibSBML(pro.setSpecies(str(product)+'__64__'+str(self.compartmentId)), 'assign product species')
            #TODO: check to see the consequences of heterologous parameters not being constant
            self._checklibSBML(pro.setConstant(True), 'set "constant" on species '+str(product))
            self._checklibSBML(pro.setStoichiometry(float(step['right'][product])),
                'set the stoichiometry ('+str(float(step['right'][product]))+')')
        #annotation
        annotation = '''<annotation>
  <rdf:RDF
  xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
  xmlns:bqbiol="http://biomodels.net/biology-qualifiers/"
  xmlns:bqmodel="http://biomodels.net/model-qualifiers/">'''
        # if the name of the species is MNX then we annotate it using MIRIAM compliance
        #TODO: need to add all known xref from different databases (not just MetaNetX)
        ############################ MIRIAM ############################
        annotation += '''
    <rdf:Description rdf:about="#'''+str(metaID or '')+'''">
      <bqbiol:is>
        <rdf:Bag>'''
        id_ident = {'mnx': 'metanetx.reaction/', 'rhea': 'rhea/', 'reactome': 'reactome/', 'bigg': 'bigg.reaction/', 'sabiork': 'sabiork.reaction/', 'ec': 'ec-code/', 'biocyc': 'biocyc/'}
        for dbId in reacXref:
            for cid in reacXref[dbId]:
                try:
                    annotation += '''
      <rdf:li rdf:resource="http://identifiers.org/'''+str(id_ident[dbId])+str(cid)+'''"/>'''
                except KeyError:
                    continue
        for ec in ecs:
            annotation += '''
      <rdf:li rdf:resource="http://identifiers.org/ec-code/'''+str(ec)+'''"/>'''
        ############################## IBISBA #########################
        #return the EC number associated with the original reaction 
        #print(step)
        annotation += '''
        </rdf:Bag>
      </bqbiol:is>
    </rdf:Description>
    <rdf:Ibisba rdf:about="#'''+str(metaID or '')+'''">
      <ibisba:ibisba xmlns:ibisba="http://ibisba.eu">
        <ibisba:smiles>'''+str(reaction_smiles or '')+'''</ibisba:smiles>
        <ibisba:rule_id>'''+str(step['rule_id'] or '')+'''</ibisba:rule_id>
        <ibisba:rule_mnxr>'''+str(step['mnxr'] or '')+'''</ibisba:rule_mnxr>
        <ibisba:rule_score value="'''+str(step['rule_score'] or '')+'''" />
        <ibisba:path_id value="'''+str(step['path_id'])+'''"/>
        <ibisba:step_id value="'''+str(step['step'])+'''"/>
        <ibisba:sub_step_id value="'''+str(step['sub_step'])+'''"/>
      </ibisba:ibisba>
    </rdf:Ibisba>
  </rdf:RDF>
</annotation>'''
        self._checklibSBML(reac.setAnnotation(annotation), 'setting annotation for reaction '+str(reacId))
        #### GROUPS #####
        if not hetero_group==None:
            newM = hetero_group.createMember()
            self._checklibSBML(newM, 'Creating a new groups member')
            self._checklibSBML(newM.setIdRef(reacId), 'Setting name to the groups member')
        elif not self.hetero_group==None:
            newM = self.hetero_group.createMember()
            self._checklibSBML(newM, 'Creating a new groups member')
            self._checklibSBML(newM.setIdRef(reacId), 'Setting name to the groups member')
        else:
            self.logger.warning('This pathway is not added to a particular group')


    ## Create libSBML reaction
    #
    # Create a reaction. fluxBounds is a list of libSBML.UnitDefinition, length of exactly 2 with the first position that is the upper bound and the second is the lower bound. reactants_dict and reactants_dict are dictionnaries that hold the following parameters: name, compartmentId, stoichiometry
    #
    # @param chemIdDictionnary containing all the cross references that we know of, can be empty)
    # @param chemXref Dictionnary containing all the cross references that we know of, can be empty
    # @param metaID Name for the reaction
    # @param inchi String Inchi associated with this species
    # @param smiles String SMILES associated with this species
    # @param compartmentId String Set this species to belong to another compartmentId than the one globally set by self.compartmentId
    # @param charge Optional parameter describing the charge of the molecule of interest
    # @param chemForm Optional chemical formulae of the substrate (not SMILES or InChI)
    # @param dG Optinal Thermodynamics constant for this species
    # @param dG_uncert Optional Uncertainty associated with the thermodynamics of the reaction 
    def createSpecies(self, 
            chemId,
            compartmentId,
            metaName=None,
            chemXref={},
            inchi=None,
            inchiKey=None,
            smiles=None,
            metaID=None):
            #TODO: add these at some point -- not very important
            #charge=0,
            #chemForm=''):
        spe = self.model.createSpecies()
        self._checklibSBML(spe, 'create species')
        ##### FBC #####
        spe_fbc = spe.getPlugin('fbc')
        self._checklibSBML(spe_fbc, 'creating this species as an instance of FBC')
        #spe_fbc.setCharge(charge) #### These are not required for FBA 
        #spe_fbc.setChemicalFormula(chemForm) #### These are not required for FBA
        #if compartmentId:
        self._checklibSBML(spe.setCompartment(compartmentId), 'set species spe compartment')
        #else:
        #    #removing this could lead to errors with xref
        #    self._checklibSBML(spe.setCompartment(self.compartmentId), 'set species spe compartment')
        #ID same structure as cobrapy
        #TODO: determine if this is always the case or it will change
        self._checklibSBML(spe.setHasOnlySubstanceUnits(False), 'set substance units')
        self._checklibSBML(spe.setBoundaryCondition(False), 'set boundary conditions')
        self._checklibSBML(spe.setConstant(False), 'set constant')
        #useless for FBA (usefull for ODE) but makes Copasi stop complaining
        self._checklibSBML(spe.setInitialConcentration(1.0), 'set an initial concentration')
        #same writting convention as COBRApy
        self._checklibSBML(spe.setId(str(chemId)+'__64__'+str(compartmentId)), 'set species id')
        if metaID==None:
            metaID = self._genMetaID(chemId)
        self._checklibSBML(spe.setMetaId(metaID), 'setting reaction metaID')

        if not metaName==None:
            self._checklibSBML(spe.setName(chemId), 'setting name for the namebolites')
        elif not metaID==None:
            self._checklibSBML(spe.setName(metaID), 'setting name for the namebolites')
        else:
            self.logger.warning('There are no inputs for the name')
        #this is setting MNX id as the name
        #this is setting the name as the input name
        
        ###### annotation ###
        annotation = '''<annotation>
  <rdf:RDF 
  xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" 
  xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" 
  xmlns:bqmodel="http://biomodels.net/model-qualifiers/">'''
        # if the name of the species is MNX then we annotate it using MIRIAM compliance
        #TODO: need to add all known xref from different databases (not just MetaNetX)
        annotation += '''
    <rdf:Description rdf:about="#'''+str(metaID or '')+'''">
      <bqbiol:is>
        <rdf:Bag>'''
        id_ident = {'mnx': 'metanetx.chemical/', 'chebi': 'chebi/CHEBI:', 'bigg': 'bigg.metabolite/', 'hmdb': 'hmdb/', 'kegg_c': 'kegg.compound/', 'kegg_d': 'kegg.drug/', 'biocyc': 'biocyc/META:', 'seed': 'seed.compound/', 'metacyc': 'metacyc/', 'sabiork': 'seed.compound/', 'reactome': 'reactome.compound/'}
        for dbId in chemXref:
            for cid in chemXref[dbId]:
                try:
                    if dbId == 'kegg' and cid[0] == 'C':
                        annotation += '''        
      <rdf:li rdf:resource="http://identifiers.org/'''+id_ident['kegg_c']+str(cid)+'''"/>'''
                    elif dbId == 'kegg' and cid[0] == 'D':
                        annotation += '''        
      <rdf:li rdf:resource="http://identifiers.org/'''+id_ident['kegg_d']+str(cid)+'''"/>'''
                    else:
                        annotation += '''        
      <rdf:li rdf:resource="http://identifiers.org/'''+str(id_ident[dbId])+str(cid)+'''"/>'''
                except KeyError:
                    continue
        annotation += '''
        </rdf:Bag>
      </bqbiol:is>
    </rdf:Description>'''   
        ###### IBISBA additional information ########
        annotation += '''
    <rdf:Ibisba rdf:about="#'''+str(metaID or '')+'''">
      <ibisba:ibisba xmlns:ibisba="http://ibisba.eu/qualifiers">
        <ibisba:smiles>'''+str(smiles or '')+'''</ibisba:smiles>
        <ibisba:inchi>'''+str(inchi or '')+'''</ibisba:inchi>
        <ibisba:inchikey>'''+str(inchiKey or '')+'''</ibisba:inchikey>
      </ibisba:ibisba>
    </rdf:Ibisba>'''
        annotation += '''
  </rdf:RDF>
</annotation>'''
        self._checklibSBML(spe.setAnnotation(annotation), 'setting the annotation for new species')

    
    ## Create libSBML pathway
    #
    # Create the collection of reactions that constitute the pathway using the Groups package and create the custom IBIBSA annotations
    #
    # @param model libSBML model to add the unit definition
    # @param reaction_id Reaction ID
    # @param name Name for the reaction
    # @param fluxBounds list of size 2 that describe the FBC upper and lower bounds for this reactions flux
    # @param reactants list of species that are the reactants of this reaction
    # @param products list of species that are the products of this reaction
    # @param reaction_smiles String smiles description of this reaction (added in IBISBA annotation)
    # @return hetero_group The number libSBML groups object to pass to createReaction to categorise the new reactions
    def createPathway(self, path_id, metaID=None, pathId='rp_pathway'):
        groups_plugin = self.model.getPlugin('groups')
        self.hetero_group = groups_plugin.createGroup()
        self.hetero_group.setId(pathId)
        if metaID==None:
            metaID = self._genMetaID(pathId)
        self.hetero_group.setMetaId(metaID)
        self.hetero_group.setKind(libsbml.GROUP_KIND_COLLECTION)
        annotation = '''<annotation>
  <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" 
  xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
    <rdf:Ibisba rdf:about="#'''+str(metaID or '')+'''">
      <ibisba:ibisba xmlns:ibisba="http://ibisba.eu">
      </ibisba:ibisba>
    </rdf:Ibisba>
  </rdf:RDF>
</annotation>'''
        self.hetero_group.setAnnotation(annotation)


    ## Create libSBML gene
    #
    # Create the list of genes in the model including its custom IBISBA annotatons
    #
    # @param model libSBML model to add the unit definition
    # @param reac libSBML reaction object
    # @param step_id The step for the number of 
    # @return libSBML gene object
    def createGene(self, reac, step_id, metaID=None):
        #TODO: pass this function to Pablo for him to fill with parameters that are appropriate for his needs
        geneName = 'RP'+str(step_id)+'_gene'
        fbc_plugin = self.model.getPlugin('fbc')
        #fbc_plugin = reac.getPlugin("fbc")
        gp = fbc_plugin.createGeneProduct()
        gp.setId(geneName)
        if metaID==None:
            metaID = self._genMetaID(str(geneName))
        gp.setMetaId(metaID)
        gp.setLabel('gene_'+str(step_id))
        gp.setAssociatedSpecies('RP'+str(step_id))
        ##### NOTE: The parameters here require the input from Pablo to determine what he needs
        annotation = '''<annotation>
  <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" 
        xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
    <rdf:Ibisba rdf:about="#'''+str(metaID or '')+'''">
      <ibisba:ibisba xmlns:ibisba="http://ibisba.eu">
        <ibisba:fasta value="" />
      </ibisba:ibisba>
    </rdf:Ibisba>
  </rdf:RDF>
</annotation>'''
        gp.setAnnotation(annotation)
    

    ## Create libSBML flux objective 
    #
    # Using the FBC package one can add the FBA flux objective directly to the model. This function sets a particular reaction as objective with maximization or minimization objectives
    #
    # @param model libSBML model to add the unit definition
    # @param fluxObjID The id given to this particular objective
    # @param reactionName The name or id of the reaction that we are setting a flux objective
    # @param coefficient FBA coefficient 
    # @param isMax Boolean to determine if we are maximizing or minimizing the objective
    # @param metaID Set the metaID
    # @return Boolean exit code
    def createFluxObj(self, fluxObjID, reactionName, coefficient, isMax=True):#, metaID=None):
        #TODO: define more complex objectives and test the simulation using a libSBML object
        fbc_plugin = self.model.getPlugin('fbc')
        target_obj = fbc_plugin.createObjective()
        target_obj.setId(fluxObjID)
        #if metaID==None:
        #    metaID = self._genMetaID(fluxObjID)
        #target_obj.setMetaId(metaID)
        if isMax:
            target_obj.setType('maximize')
        else:
            target_obj.setType('minimize')
        fbc_plugin.setActiveObjectiveId(fluxObjID) # this ensures that we are using this objective when multiple
        target_flux_obj = target_obj.createFluxObjective()
        target_flux_obj.setReaction(reactionName)
        target_flux_obj.setCoefficient(coefficient)


    ## Generate a generic model 
    #
    # Since we will be using the same type of parameters for the RetroPath model, this function
    # generates a libSBML model with parameters that will be mostly used
    #
    #
    #
    def genericModel(self, modelName, modelID, compXref, compartment_id):
        self.createModel(modelName, modelID)
        # mmol_per_gDW_per_hr
        unitDef = self.createUnitDefinition('mmol_per_gDW_per_hr')
        self.createUnit(unitDef, libsbml.UNIT_KIND_MOLE, 1, -3, 1)
        self.createUnit(unitDef, libsbml.UNIT_KIND_GRAM, 1, 0, 1)
        self.createUnit(unitDef, libsbml.UNIT_KIND_SECOND, 1, 0, 3600)
        # kj_per_mol
        gibbsDef = self.createUnitDefinition('kj_per_mol')
        self.createUnit(gibbsDef, libsbml.UNIT_KIND_JOULE, 1, 3, 1)
        self.createUnit(gibbsDef, libsbml.UNIT_KIND_MOLE, -1, 1, 1)
        # infinity parameters (FBA)
        #upInfParam = self.createParameter('B_INF', float('inf'), 'kj_per_mol')
        #lowInfParam = self.createParameter('B__INF', float('-inf'), 'kj_per_mol')
        upNineParam = self.createParameter('B__999999', -999999, 'mmol_per_gDW_per_hr')
        lowNineParam = self.createParameter('B_999999', 999999, 'mmol_per_gDW_per_hr')
        lowZeroParam = self.createParameter('B_0', 0, 'mmol_per_gDW_per_hr')
        #compartment
        #TODO: create a new compartment 
        #self.createCompartment(1, 'MNXC3', 'cytoplasm', compXref)
        #try to recover the name from the Xref
        try:
            name = compXref['name'][0]
        except KeyError:
            name = compartment_id+'_name'
        self.createCompartment(1, compartment_id, name, compXref)


    ##################################################################################################
    ########################################### TESTING ##############################################
    ##################################################################################################


    ## Test writing of an SBML model
    #
    # Private test function that will write a given pathway to file
    #
    # @param cofactors_rp_paths Dictionnary of the outputs from RP2paths
    # @param path_id Integer The cofactors_rp_paths pathway to write to SBML
    def _testWriter(self, cofactors_rp_paths=None, path_id=None):
        #### Define the if global or local parameter is to be used
        try:
            if cofactors_rp_paths==None and self.cofactors_rp_paths==None:
                raise TypeError
            if cofactors_rp_paths==None and not self.cofactors_rp_paths==None:
                cofactors_rp_paths = self.cofactors_rp_paths
        except EmptyOutRPpaths:
            self.logger.error('Both class and input cofactors_rp_paths are empty')
        #### Check what is the best path to use
        if path_id==None:
            path_id = 1
        return None

if __name__ == "__main__":
    #read the TAR.XZ with all the SBML pathways
    rpsbml_paths = {}
    tar = tarfile.open('tests/testFBAin.tar.xz') #TODO: create this
    rpsbml_paths = {}
    for member in tar.getmembers():
        rpsbml_paths[member.name] = rpFBA.rpSBML(member.name,libsbml.readSBMLFromString(tar.extractfile(member).read().decode("utf-8")))
    #pass the different models to the SBML solvers and write the results to file
    #TODO
    #designed to write using TAR.XZ with all the SBML pathways
    with tarfile.open('testFBAout.tar.xz', 'w:xz') as tf:
        for rpsbml_name in rpsbml_paths:
            data = libsbml.writeSBMLToString(rpsbml_paths[rpsbml_name].document).encode('utf-8')
            fiOut = BytesIO(data)
            info = tarfile.TarInfo(rpsbml_name)
            info.size = len(data)
            tf.addfile(tarinfo=info, fileobj=fiOut)


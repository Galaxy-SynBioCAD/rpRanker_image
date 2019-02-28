import libsbml
from hashlib import md5
import logging
import os

class Annotate():
    """ Class to add, modify or delete IBIBSA annotations and other FBA parameters
    """

    def __init__(self, path=None):
        self.path = path
        self.model = None

    ####### helper 

    def _exportSBML(self, type_name, model, model_id, path=None):
        """Export the libSBML model to an SBML file
        """
        ####### check the path #########
        p = None
        if path:
            if path[-1:]=='/':
                path = path[:-1]
            if not os.path.isdir(path):
                if outputPath:
                    p = outputPath
                else:
                    logging.error('The output path is not a directory: '+str(path))
                    return False
            else:
                p = path
        else:
            p = outputPath
        p += '/'+type_name
        ########## check and create folder #####
        if not os.path.exists(p):
            os.makedirs(p)
        else:
            logging.error('Cannot recognise input type')
            return False
        libsbml.writeSBMLToFile(model, p+'/'+str(model_id)+'.sbml')
        return True


    def _openSBML(self, ):
        """Situation where an SBML is passed to add the heterologous pathway
        """
    


    def _checklibSBML(self, value, message):
        """Check that the libSBML python calls do not return error INT and if so, display the error
        Taken from: http://sbml.org/Software/libSBML/docs/python-api/create_simple_model_8py-example.html
        """
        if value is None:
            logging.error('LibSBML returned a null value trying to ' + message + '.')
            raise SystemExit('LibSBML returned a null value trying to ' + message + '.')
        elif type(value) is int:
            if value==libsbml.LIBSBML_OPERATION_SUCCESS:
                return
            else:
                err_msg = 'Error encountered trying to ' + message + '.' \
                        + 'LibSBML returned error code ' + str(value) + ': "' \
                        + libsbml.OperationReturnValue_toString(value).strip() + '"'
                logging.error(err_msg)
                raise SystemExit(err_msg)
        else:
            logging.info(message)
            return

    def _nameToSbmlId(self, name):
        """Function to rewrite a string to libSBML metaid valid format
        """
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


    def _genMetaID(self, name):
        """Generate the metaID from the name of the parameter
        """
        return nameToSbmlId(md5(str(name).encode('utf-8')).hexdigest())        


    def createModel(self, name, model_id):
        """Function that creates a new libSBML model instance and initiates it with the appropriate
        packages. Creates a cytosol compartment
        """
        ## sbmldoc
        sbmlns = libsbml.SBMLNamespaces(3,1)
        sbmlns.addPkgNamespace('groups',1)
        sbmlns.addPkgNamespace('fbc',2)
        #sbmlns = libsbml.SBMLNamespaces(3,1,'groups',1)
        sbmlDoc = libsbml.SBMLDocument(sbmlns)
        sbmlDoc.setPackageRequired('fbc', False) #!!!! must be set to false for no apparent reason
        sbmlDoc.setPackageRequired('groups', False) #!!!! must be set to false for no apparent reason
        ## sbml model
        model = sbmlDoc.createModel()
        model.setId(model_id)
        model.setName(name)
        checklibSBML(self.model, 'create model')
        checklibSBML(self.model.setTimeUnits('second'), 'setting model time unit')
        checklibSBML(self.model.setExtentUnits('mole'), 'setting model compartment unit')
        checklibSBML(self.model.setSubstanceUnits('mole'), 'setting model substance unit')
        return True

    def createCompartment(self, model, size, name, compartment_id):
        ## cytoplasm compartment TODO: consider seperating it in another function 
        # if another compartment is to be created
        comp = model.createCompartment()
        checklibSBML(comp, 'create compartment')
        checklibSBML(comp.setId(compartment_id), 'set compartment id')
        checklibSBML(comp.setConstant(True), 'set compartment "constant"')
        checklibSBML(comp.setSize(size), 'set compartment "size"')
        checklibSBML(comp.setSBOTerm(290), 'set SBO term for the cytoplasm compartment')
        checklibSBML(comp.setName(name), 'set the name for the cytoplam')
        return True

    def createUnitDefinition(self, model, unit_id):
        """Function that creates a unit definition (composed of one or more units)
        """
        unitDef = model.createUnitDefinition()
        checklibSBML(unitDef, 'creating flux unit definition')
        checklibSBML(unitDef.setId(unit_id), 'setting flux id')
        checklibSBML(unitDef.setMetaId(genMetaID(unit_id)), 'Setting flux metaID')
        return True

    def createUnit(self, unitDef, libsbmlunit, exponent, scale, multiplier):
        """Function that created a unit
        """
        unit = unitDef.createUnit()
        checklibSBML(unit, 'creating unit')
        checklibSBML(unit.setKind(libsbmlunit), 'setting the kind of mole unit')
        checklibSBML(unit.setExponent(exponent), 'setting the exponenent of the mole unit')
        checklibSBML(unit.setScale(scale), 'setting the scale of the mole unit')
        checklibSBML(unit.setMultiplier(multiplier), 'setting the multiplier of the mole unit')
        return True

    def createParameter(self, parameter_id, value, unit):
        """Parameters, in our case, used for the bounds for FBA analysis.
        unit parameter must be an instance of unitDefinition
        """
        infParam = self.model.createParameter()
        checklibSBML(infParam, 'Creating a new parameter object')
        checklibSBML(infParam.setConstant(True), 'setting INF as constant')
        checklibSBML(infParam.setId(parameter_id), 'setting INF ID')
        checklibSBML(infParam.setValue(value), 'setting value of INF')
        checklibSBML(infParam.setUnits(unit), 'setting INF units')
        checklibSBML(infParam.setSBOTerm(625), 'setting INF SBO term')
        checklibSBML(infParam.setMetaId(genMetaID(parameter_id), 'setting INF meta ID')
        return True

    def createReaction(self, model, reaction_id, name, fluxBounds, reactants, products, reaction_smiles):
        """Create a reaction. fluxBounds is a list of libSBML.UnitDefinition, length of exactly 2 with the 
        first position that is the upper bound and the second is the lower bound. reactants_dict and 
        reactants_dict are dictionnaries that hold the following parameters: name, compartments, stoichiometry 
        """
        reac = model.createReaction()
        checklibSBML(reac, 'create reaction')
        ################ FBC ####################
        reac_fbc = reac.getPlugin('fbc')
        checklibSBML(reac_fbc, 'extending reaction for FBC')
        #bounds
        if len(fluxBunds)==2 and isinstance(fluxBounds, list):
            checklibSBML(reac_fbc.setUpperFluxBound(fluxBounds[0]), 'setting '+name+' upper flux bound')
            checklibSBML(reac_fbc.setLowerFluxBound(fluxBounds[1]), 'setting '+name+' lower flux bound')
        else:
            logging.error('fluxBounds is either not len()==2 or is not a list')
            return False
        #########################################
        #reactions
        checklibSBML(reac.setId(reaction_id), 'set reaction id') #same convention as cobrapy
        checklibSBML(reac.setName(name), 'set name') #same convention as cobrapy
        checklibSBML(reac.setSBOTerm(185), 'setting the system biology ontology (SBO)') #set as process
        checklibSBML(reac.setReversible(True), 'set reaction reversibility flag')
        checklibSBML(reac.setFast(False), 'set reaction "fast" attribute')
        metaID = genMetaID(reaction_id)
        checklibSBML(reac.setMetaId(metaID), 'setting species metaID')
        #reactants_dict
        for r in reactants:
            spe = reac.createReactant()
            checklibSBML(spe, 'create reactant')
            #use the same writing convention as CobraPy
            checklibSBML(spe.setSpecies(str(r['name'])+'__64__'+str(r['compartment'])), 'assign reactant species')
            #TODO: check to see if this must not be an input
            checklibSBML(spe.setConstant(True), 'set "constant" on species '+str(r['name']))
            checklibSBML(spe.setStoichiometry(abs(float(r['stoichiometry']))), 
                'set stoichiometry ('+str(abs(float(r['stoichiometry'])))+')')
        #reactants_dict
        for p in products:
            pro = reac.createProduct()
            checklibSBML(pro, 'create product '+str(meta))
            checklibSBML(pro_r.setSpecies(str(meta)+'__64__'+str(compartment)), 'assign product species')
            #check to see that this does not need to be an input
            checklibSBML(pro_r.setConstant(True), 'set "constant" on species '+str(meta))
            checklibSBML(pro_r.setStoichiometry(float(p['stoichiometry'])), 
                'set the stoichiometry ('+str(float(p['stoichiometry']))+')')
        #annotation
        annotation = '''<annotation>
  <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
  xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
    <rdf:Description rdf:about="#'''+str(metaID)+'''">
      <ibisba:ibisba xmlns:ibisba="http://ibisba.eu">
        <ibisba:smiles>'''+str(reaction_smiles)+'''</ibisba:smiles>
        <ibisba:parameter type="ddG" units="kj_per_mol" value="None"/>
        <ibisba:parameter type="ddG_uncert" units="kj_per_mol" value="None"/>
      </ibisba:ibisba>
    </rdf:Description>
  </rdf:RDF>
</annotation>'''
        checklibSBML(reac.setAnnotation(annotation), 'setting annotation for reaction '+str(name))
        return True


    def createSpecies(self, model, compartment, name):
        """Create a libSBML species with FBC parameters and custom IBISBA annotations 
        """
        spe = model.createSpecies()
        checklibSBML(spe, 'create species')
        ##### FBC #####
        spe_fbc = spe.getPlugin('fbc')
        checklibSBML(spe_fbc, 'creating this species as an instance of FBC')
        spe_fbc.setCharge(0) #### TO BE DETERMINED
        spe_fbc.setChemicalFormula('') #### TO BE DETERMINED
        checklibSBML(spe.setCompartment(compartment), 'set species spe compartment')
        #ID same structure as cobrapy
        #TODO: determine if this is always the case or it will change
        checklibSBML(spe.setHasOnlySubstanceUnits(False), 'set substance units')
        checklibSBML(spe.setBoundaryCondition(False), 'set boundary conditions')
        checklibSBML(spe.setConstant(False), 'set constant')
        #useless but makes Copasi stop complaining
        checklibSBML(spe.setInitialConcentration(1.0), 'set an initial concentration')
        #same writting convention as COBRApy
        checklibSBML(spe.setId(nameToSbmlId(str(name)+'__64__'+str(compartment))), 'set species id')
        metaID = _genMetaID(name)
        checklibSBML(spe.setMetaId(metaID), 'setting reaction metaID')
        checklibSBML(spe.setName(name), 'setting name for the metabolites')
        ###### annotation ###
        annotation = '''<annotation>
  <rdf:RDF 
  xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" 
  xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" 
  xmlns:bqmodel="http://biomodels.net/model-qualifiers/"
  xmlns:ibisba="http://ibisba.eu/qualifiers">'''
        # if the name of the species is MNX then we annotate it using MIRIAM compliance
        if meta[:3]=='MNX':
            annotation += '''
  <rdf:Description rdf:about="#'''+str(metaID)+'''">
    <bqbiol:is>
      <rdf:Bag>
        <rdf:li rdf:resource="http://identifiers.org/metanetx.chemical/'''+str(meta)+'''"/>
      </rdf:Bag>
    </bqbiol:is>
  </rdf:Description>'''   
        #add IBISBA additional information
        if meta_smiles[meta]:
            annotation += '''
  <rdf:Description rdf:about="#'''+str(metaID)+'''">
    <ibisba:ibisba xmlns:ibisba="http://ibisba.eu/qualifiers">
      <ibisba:smiles>'''+str(meta_smiles[meta])+'''</ibisba:smiles>'''
        else:
            annotation += '''
            <rdf:Description rdf:about="#'''+str(metaID)+'''">
              <ibisba:ibisba xmlns:ibisba="http://ibisba.eu/qualifiers">
                <ibisba:smiles></ibisba:smiles>'''
        if meta_inchi[meta]:
            annotation += '''
            <ibisba:inchi>'''+str(meta_inchi[meta])+'''</ibisba:inchi>
            <ibisba:inchikey>'''+str(Chem.rdinchi.InchiToInchiKey(meta_inchi[meta]))+'''</ibisba:inchikey>'''
        else:
            annotation += '''
            <ibisba:inchi></ibisba:inchi>
            <ibisba:inchikey></ibisba:inchikey>'''
        annotation += '''
            <ibisba:parameter type="ddG" units="kj_per_mol" value="None"/>
            <ibisba:parameter type="ddG_uncert" units="kj_per_mol" value="None"/>
          </ibisba:ibisba>
        </rdf:Description>'''
        annotation += '''
      </rdf:RDF>
    </annotation>'''
        checklibSBML(spe.setAnnotation(annotation), 'setting the annotation for new species')
    


    def createPathway(self, ):
        """Create the collection of reactions that constitute the pathway using the Groups
        package and create the custom IBIBSA annotations
        """
        

    def createGene(self, ):
        """Create the list of genes in the model including its custom IBISBA annotatons 
        """

    
    def createFluxObj(self, ):
    

    #TODO: write the function but seems like an overkill since it makes more sense to define the boundaries
    #at the creation of the reactions
    def createFluxBounds(self, ):

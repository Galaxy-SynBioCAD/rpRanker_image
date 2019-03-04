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


    ##########################################
    ####### Private functions ################ 
    ##########################################


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


    def _openSBML(self, filename):
        """Situation where an SBML is passed to add the heterologous pathway
        """
        document = libSBML.readSBML(filename)
        errors = document.getNumErrors()
        if errors>0:
            logger.warning('Reading the document has returned some errors ('+str(errors)+')')
            return False
        return document.model


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


    ######################################
    ############# WRITE ##################
    ######################################


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


    def createReaction(self, model, reaction_id, name, reactants, products, reaction_smiles):
        """Create a reaction. fluxBounds is a list of libSBML.UnitDefinition, length of exactly 2 with the 
        first position that is the upper bound and the second is the lower bound. reactants_dict and 
        reactants_dict are dictionnaries that hold the following parameters: name, compartments, stoichiometry 
        """
        reac = model.createReaction()
        checklibSBML(reac, 'create reaction')
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
        <ibisba:ddG units="kj_per_mol" value="None"/>
        <ibisba:ddG_uncert units="kj_per_mol" value="None"/>
      </ibisba:ibisba>
    </rdf:Description>
  </rdf:RDF>
</annotation>'''
        checklibSBML(reac.setAnnotation(annotation), 'setting annotation for reaction '+str(name))
        return True


    def createSpecies(self, model, compartment, name, group=None, smiles=None, inchi=None, inchikey=None):
        """Create a libSBML species with FBC parameters and custom IBISBA annotations. Note that this
        assumes that the species is newly created
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
        if name[:3]=='MNX':
            annotation += '''
    <rdf:Description rdf:about="#'''+str(metaID)+'''">
      <bqbiol:is>
        <rdf:Bag>
          <rdf:li rdf:resource="http://identifiers.org/metanetx.chemical/'''+str(name)+'''"/>
        </rdf:Bag>
      </bqbiol:is>
    </rdf:Description>'''   
        #add IBISBA additional information
        ##########SMILES#############
        if smiles:
            annotation += '''
    <rdf:Description rdf:about="#'''+str(metaID)+'''">
      <ibisba:ibisba xmlns:ibisba="http://ibisba.eu/qualifiers">
        <ibisba:smiles> value="'''+str(smiles)+'''" />'''
        else:
            annotation += '''
    <rdf:Description rdf:about="#'''+str(metaID)+'''">
      <ibisba:ibisba xmlns:ibisba="http://ibisba.eu/qualifiers">
        <ibisba:smiles value="" />'''
        ################INCHI############
        if inchi:
            annotation += '''
        <ibisba:inchi> value="'''+str(inchi)+'''" />'''
        else:
            annotation += '''
        <ibisba:inchi value="" />'''
        ###############INCHIKEY###########
        if inchikey:
            annotation += '''
        <ibisba:inchikey value="'''+str(inchikey)+'''" />'''
        elif inchi and not inchikey:
            annotation += '''
        <ibisba:inchikey value="'''+str(Chem.rdinchi.InchiToInchiKey(meta_inchi[meta]))+'''" />'''
        else:
            annotation += '''
        <ibisba:inchikey value="" />'''
        annotation += '''
        <ibisba:ddG units="kj_per_mol" value="None"/>
        <ibisba:ddG_uncert units="kj_per_mol" value="None"/>
      </ibisba:ibisba>
    </rdf:Description>'''
        annotation += '''
  </rdf:RDF>
</annotation>'''
        checklibSBML(spe.setAnnotation(annotation), 'setting the annotation for new species')
        return True    


    def createPathway(self, model, pathway_id, name, reaction_refIDs):
        """Create the collection of reactions that constitute the pathway using the Groups
        package and create the custom IBIBSA annotations
        The metaID_list is a list of strings corresponding to the metaID of the reactions that
        are to be added to the 
        """
        groups_plugin = model.getPlugin("groups")
        hetero_group = groups_plugin.createGroup()
        checklibSBML(hetero_group, 'creating a groups parameter')
        checklibSBML(hetero_group.setId(pathway_id), 'setting pathway id')
        metaID = _genMetaID(pathway_id)
        checklibSBML(hetero_group.setMetaId(metaID), 'setting pathway metaID')
        checklibSBML(hetero_group.setName(name), 'setting pathway name')
        #we should always need to define this group as a COLLECTION
        checklibSBML(hetero_group.setKind(libsbml.GROUP_KIND_COLLECTION), 'setting the groups as a collection')
        annotation = '''<annotation>
  <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" 
  xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
    <rdf:Description rdf:about="#'''+str(metaID)+'''">
      <ibisba:ibisba xmlns:ibisba="http://ibisba.eu">
        <ibisba:ddG units="kj_per_mol" value="None"/>
        <ibisba:ddG_uncert units="kj_per_mol" value="None"/>
      </ibisba:ibisba>
    </rdf:Description>
  </rdf:RDF>
</annotation>'''
        hetero_group.setAnnotation(annotation)
        for r in reaction_refIDs:
            newMember = group.createMember()
            checklibSBML(newMember, 'creating a new groups member')
            checklibSBML(newMember.setIdRef(r))
        return True


    def createGene(self, model, name, gene_id, label, associated_reaction):
        """Create the list of genes in the model including its custom IBISBA annotatons. Note that
        this assumes the creation of a gene that does not exist and will overwrite an existing one
        """
        fbc_plugin = model.getPlugin('fbc')
        gp = fbc_plugin.createGeneProduct()
        checklibSBML(gp, 'creating gene product')
        checklibSBML(gp.setId(name), 'setting gene name')
        metaID = _genMetaID(gene_id)
        checklibSBML(gp.setMetaId(metaID), 'setting gene metaID')
        checklibSBML(gp.setLabel(label), 'setting gene label')
        checklibSBML(gp.setAssociatedSpecies(associated_reaction), 'setting gene associated reaction')
        annotation = '''<annotation>
  <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" 
  xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
    <rdf:Description rdf:about="#'''+str(metaID)+'''">
      <ibisba:ibisba xmlns:ibisba="http://ibisba.eu">
        <ibisba:fasta value="" />
      </ibisba:ibisba>
    </rdf:Description>
  </rdf:RDF>
</annotation>'''
        checklibSBML(gp.setAnnotation(annotation), 'setting gene annotation')
        return True


    #note that this might need to be on a case by case basis, since we have the ability
    #to define multiple objectives
    def createFluxObj(self, model, name, fluxObj_id, obj_type, associated_reaction)
        """Create the flux objective for the model
        """
        fbc_plugin = model.getPlugin("fbc")
        target_obj = fbc_plugin.createObjective()
        target_obj.setId(fluxObj_id)
        target_obj.setMetaId(_genMetaID(fluxObj_id))
        if obj_type=='maximize' or obj_type=='minimize':
            target_obj.setType(obj_type)
        else:
            logger.error('the input flux objective must be either maximize or minimize')
            return False
        fbc_plugin.setActiveObjectiveId('target_obj') # this ensures that we are using this objective when multiple
        target_flux_obj = target_obj.createFluxObjective()
        target_flux_obj.setReaction(associated_reaction)
        target_flux_obj.setCoefficient(1)
        return True


    #TODO: write the function but seems like an overkill since it makes more sense to define the boundaries
    #at the creation of the reactions
    def createFluxBounds(self, model, name, reaction_id, upperBound, lowerBound):
        """Create the flux bounds for an input reaction (reaction_id)
        """
        #TODO: wrap it with try/catch in case the reaction does not exist
        reac = model.getReaction(reaction_id)
        reac_fbc = reac.getPlugin('fbc')
        checklibSBML(reac_fbc, 'extending reaction for FBC')
        #bounds
        checklibSBML(reac_fbc.setUpperFluxBound(upperBound), 'setting '+name+' upper flux bound')
        checklibSBML(reac_fbc.setLowerFluxBound(lowerBound), 'setting '+name+' lower flux bound')
        return True 

    #########################################
    ################ READ ###################
    #########################################


    #NOTE: these are helper functions, and the native libsbml functions should be used preferably
    
    def readAnnotation(self, annotation):
        """Read the annotations of a reaction or a species etc... Structure should be the same, for each. That is:
        level 1: annotation
        level 2: rdf:RDF
        level 3: rdf:Description
        level 4: bqbiol:is/ibisba:ibisba --> former is the MIRIAM species and the other are the 
        nest levels are the same
        """
        annotation = sbase.getAnnotation()
        toRet_annot = {'miriam': {}, 'ibisba': {}}
        if annotation.hasChild('RDF') and annotation.getChild('RDF').hasChild('Description'):
            #test that there is at least one
            #must consider the case that the first IBIBSA and MIRIAM annotations are "swapped"
            for i in annotation.getChild('RDF').getNumChildren():
                if annotation.getChild('RDF').getChild(i).hasChild('is'):
                    for i in annotation.getChild('RDF').getChild(0).getChild('is').getChild('Bag').getNumChildren():
                        m_a = toRet_annot.getChild('RDF').getChild(0).getChild('is').getChild('Bag').getChild(i).getAttrValue().split('/')
                        #should we be splitting this here?
                        if not toRet_annot['miriam'][m_a[-2].split('.')[0]]:
                            toRet_annot['miriam'][m_a[-2].split('.')[0]] = []
                        toRet_annot['miriam'][m_a[-2].split('.')[0]].append(m_a[-1].splt(':')[1])
                elif annotation.getChild('RDF').getChild(i).hasChild('ibibsa'):
                    i_a = toRet_annot.getChild('RDF').getChild(0).getChild('ibisba')
                    #WARNING: cannot check that they exist using the 
                    #TODO: wrap this around try/catch in case some don't have some of these
                    toRet_annot['ibisba']['smiles'] = i_a.getChild('smiles').getAttrValue('value')
                    toRet_annot['ibisba']['inchi'] = i_a.getChild('inchi').getAttrValue('value')
                    toRet_annot['ibisba']['inchikey'] = i_a.getChild('inchikey').getAttrValue('value')
                    toRet_annot['ibisba']['ddG'] = i_a.getChild('ddG').getAttrValue('value')
                    toRet_annot['ibisba']['ddG_uncert'] = i_a.getChild('ddG_uncert').getAttrValue('value')
        else:
            logging.error('Either the structure is wrong or the annoation for the passed SBase is empty')



    #Need to define the datatype that can be passed between the different reader functions
    def readGenes(self, model):
        toRet = []
        
            gene = {'metaid': , 'fbc_id': , 'fbc_label': , 'fbc_type': , 'flux_obj': []}


    def readPathway(self, model):
        """Helper function to read the pathway. Note that these are not required and one is advised to use
        native libsbml functions to read the required information directly from the model object
        """
        groups_plugin = model.getPlugin('groups')
        rp_pathway = groups_plugin.getGroup('rp_pathway')
        annotation = rp_pathway.getAnnotation()
        return readAnnotation(annotation)['ibisba']




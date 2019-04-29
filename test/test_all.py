import sys
import os
import libsbml
import pytest
import pickle
import gzip

sys.path.append(os.path.join(os.path.dirname('__file__'), '..'))
import rpSBML
#import ../rpSBML


path_id = 1
steps = [{'right': {'CMPD_0000000003': 1, 'MNXM13': 1, 'MNXM15': 1, 'MNXM8': 1}, 'left': {'MNXM10': 1, 'MNXM188': 1, 'MNXM4': 1, 'MNXM1': 3}}, {'right': {'TARGET_0000000001': 1, 'MNXM1': 2}, 'left': {'CMPD_0000000003': 1, 'MNXM4': 1, }}]
reaction_smiles = ['[H]Oc1c([H])c([H])c([H])c([H])c1O[H]>>O=O.[H]N=C(O[H])C1=C([H])N(C2([H])OC([H])(C([H])([H])OP(=O)(O[H])OP(=O)(O[H])OC([H])([H])C3([H])OC([H])(n4c([H])nc5c(N([H])[H])nc([H])nc54)C([H])(O[H])C3([H])O[H])C([H])(O[H])C2([H])O[H])C([H])=C([H])C1([H])[H].[H]OC(=O)c1c([H])c([H])c([H])c([H])c1N([H])[H]', '[H]OC(=O)C([H])=C([H])C([H])=C([H])C(=O)O[H]>>O=O.[H]Oc1c([H])c([H])c([H])c([H])c1O[H]', ]
rp_smiles = {'MNXM10': '[H]N=C(O[H])C1=C([H])N(C2([H])OC([H])(C([H])([H])OP(=O)(O[H])OP(=O)(O[H])OC([H])([H])C3([H])OC([H])(n4c([H])nc5c(N([H])[H])nc([H])nc54)C([H])(O[H])C3([H])O[H])C([H])(O[H])C2([H])O[H])C([H])=C([H])C1([H])[H]',
'MNXM188': '[H]OC(=O)c1c([H])c([H])c([H])c([H])c1N([H])[H]',
'MNXM4': 'O=O',
'MNXM1': '[H+]',
'CMPD_0000000003': '[H]Oc1c([H])c([H])c([H])c([H])c1O[H]',
'MNXM13': 'O=C=O',
'MNXM15': '[H]N([H])[H]',
'MNXM8': 'NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H](COP([O-])(=O)OP([O-])(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O',
'TARGET_0000000001': '[H]OC(=O)C([H])=C([H])C([H])=C([H])C(=O)O[H]'}
rp_inchi = {'MNXM10': 'InChI=1S/C21H29N7O14P2/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(32)14(30)11(41-21)6-39-44(36,37)42-43(34,35)38-5-10-13(29)15(31)20(40-10)27-3-1-2-9(4-27)18(23)33/h1,3-4,7-8,10-11,13-16,20-21,29-32H,2,5-6H2,(H2,23,33)(H,34,35)(H,36,37)(H2,22,24,25)/p-2/t10-,11-,13-,14-,15-,16-,20-,21-/m1/s1',
'MNXM188': 'InChI=1S/C7H7NO2/c8-6-4-2-1-3-5(6)7(9)10/h1-4H,8H2,(H,9,10)/p-1',
'MNXM4': 'InChI=1S/O2/c1-2',
'MNXM1': 'InChI=1S',
'CMPD_0000000003': None,
'MNXM13': 'InChI=1S/CO2/c2-1-3',
'MNXM15': 'InChI=1S/H3N/h1H3/p+1',
'MNXM8': 'InChI=1S/C21H27N7O14P2/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(32)14(30)11(41-21)6-39-44(36,37)42-43(34,35)38-5-10-13(29)15(31)20(40-10)27-3-1-2-9(4-27)18(23)33/h1-4,7-8,10-11,13-16,20-21,29-32H,5-6H2,(H5-,22,23,24,25,33,34,35,36,37)/p-1/t10-,11-,13-,14-,15-,16-,20-,21-/m1/s1',
'TARGET_0000000001': 'InChI=1S/C6H6O4/c7-5(8)3-1-2-4-6(9)10/h1-4H,(H,7,8)(H,9,10)'}

compXref = pickle.load(gzip.open(os.path.join(os.path.abspath('..'), 'cache/compXref.pickle.gz'), 'rb'))
chemXref = pickle.load(gzip.open(os.path.join(os.path.abspath('..'), 'cache/chemXref.pickle.gz'), 'rb'))
reacXref = pickle.load(gzip.open(os.path.join(os.path.abspath('..'), 'cache/reacXref.pickle.gz'), 'rb'))


class TestClass(object):

    #function to compare two models 
    #TODO: need to improve this function to check in a lower level that two SBML are the same
    # or use another strategy or package
    def compareModels(self, modelOne, modelTwo):
        try:
            for one_species in modelOne.getListOfSpecies():
                if modelTwo.getSpecies(one_species.getId()) is None:
                    return False
        except AttributeError:
            pass
        try:
            for one_reaction in modelOne.getListOfReactions():
                if modelTwo.getReaction(one_reaction.getId()) is None:
                    return False
        except AttributeError:
            pass
        try:
            for one_compartment in modelOne.getListOfCompartments():
                if modelTwo.getCompartment(one_compartment.getId()) is None:
                    return False
        except AttributeError:
            pass
        try:
            for one_parameter in modelOne.getListOfParameters():
                if modelTwo.getParameter(one_parameter.getId()) is None:
                    return False
        except AttributeError:
            pass
        try:
            for one_unitDef in modelOne.getListOfUnitDefinitions():
                if modelTwo.getUnitDefinition(one_unitDef.getId()) is None:
                    return False
        except AttributeError:
            pass
        return True



    def test_createModel(self):
        rpsbml = rpSBML.rpSBML(None, None, None, None, os.path.abspath('../cache'))
        rpsbml.createModel('RetroPath_heterologous_pathway', 'rp_model')
        inModel = libsbml.readSBML(os.path.join(os.path.abspath('test_models'), 'test_createModel.sbml'))
        assert self.compareModels(inModel.getModel(), rpsbml.model)==True


    def test_genericModel(self):
        rpsbml = rpSBML.rpSBML(None, None, None, None, os.path.abspath('../cache'))
        rpsbml.genericModel('RetroPath_heterologous_pathway', 'rp_model', compXref)
        inModel = libsbml.readSBML(os.path.join(os.path.abspath('test_models'), 'test_genericModel.sbml'))
        assert self.compareModels(inModel.getModel(), rpsbml.model)==True


    def test_createUnitDefinition(self):
        rpsbml = rpSBML.rpSBML(None, None, None, None, os.path.abspath('../cache'))
        rpsbml.createModel('RetroPath_heterologous_pathway', 'rp_model')
        unitDef = rpsbml.createUnitDefinition('mmol_per_gDW_per_hr')
        inModel = libsbml.readSBML(os.path.join(os.path.abspath('test_models'), 'test_createUnitDefinition.sbml'))
        assert self.compareModels(inModel.getModel(), rpsbml.model)==True


    def test_createUnit(self):
        rpsbml = rpSBML.rpSBML(None, None, None, None, os.path.abspath('../cache'))
        rpsbml.createModel('RetroPath_heterologous_pathway', 'rp_model')
        unitDef = rpsbml.createUnitDefinition('mmol_per_gDW_per_hr')
        rpsbml.createUnit(unitDef, libsbml.UNIT_KIND_MOLE, 1, -3, 1)
        rpsbml.createUnit(unitDef, libsbml.UNIT_KIND_GRAM, 1, 0, 1)
        rpsbml.createUnit(unitDef, libsbml.UNIT_KIND_SECOND, 1, 0, 3600)
        inModel = libsbml.readSBML(os.path.join(os.path.abspath('test_models'), 'test_createUnit.sbml'))
        assert self.compareModels(inModel.getModel(), rpsbml.model)==True


    def test_createParameter(self):
        rpsbml = rpSBML.rpSBML(None, None, None, None, os.path.abspath('../cache'))
        rpsbml.createModel('RetroPath_heterologous_pathway', 'rp_model')
        upInfParam = rpsbml.createParameter('B_999999', 999999.0, 'kj_per_mol')
        inModel = libsbml.readSBML(os.path.join(os.path.abspath('test_models'), 'test_createParameter.sbml'))
        assert self.compareModels(inModel.getModel(), rpsbml.model)==True


    def test_createCompartent(self):
        rpsbml = rpSBML.rpSBML(None, None, None, None, os.path.abspath('../cache'))
        rpsbml.createModel('RetroPath_heterologous_pathway', 'rp_model')
        rpsbml.createCompartment(1, 'MNXC3', 'cytoplasm', compXref)
        inModel = libsbml.readSBML(os.path.join(os.path.abspath('test_models'), 'test_createCompartment.sbml'))
        assert self.compareModels(inModel.getModel(), rpsbml.model)==True


    def test_createSpecies(self):
        rpsbml = rpSBML.rpSBML(None, None, None, None, os.path.abspath('../cache'))
        rpsbml.genericModel('RetroPath_heterologous_pathway', 'rp_model', compXref)
        for meta in set([i for step in steps for lr in ['left', 'right'] for i in step[lr]]):
            try:
                inchi = rp_inchi[meta]
            except KeyError:
                inchi = None
            try:
                smiles = rp_smiles[meta]
            except KeyError:
                smiles = None
            rpsbml.createSpecies(meta, chemXref, None, inchi, smiles, 'MNXC3', 0, '')
        inModel = libsbml.readSBML(os.path.join(os.path.abspath('test_models'), 'test_createSpecies.sbml'))
        assert self.compareModels(inModel.getModel(), rpsbml.model)==True


    def test_createPathway(self):
        rpsbml = rpSBML.rpSBML(None, None, None, None, os.path.abspath('../cache'))
        rpsbml.genericModel('RetroPath_heterologous_pathway', 'rp_model', compXref)
        for meta in set([i for step in steps for lr in ['left', 'right'] for i in step[lr]]):
            try:
                inchi = rp_inchi[meta]
            except KeyError:
                inchi = None
            try:
                smiles = rp_smiles[meta]
            except KeyError:
                smiles = None
            rpsbml.createSpecies(meta, chemXref, None, inchi, smiles, 'MNXC3', 0, '')
        rpsbml.createPathway(path_id)
        #reactions
        step_id = 0
        for stepNum in range(len(steps)):
            rpsbml.createReaction('RP_'+str(stepNum), 'B_999999', 'B__999999', steps[stepNum], reaction_smiles[stepNum], reacXref, None, rpsbml.hetero_group, None)
            step_id += 1
        inModel = libsbml.readSBML(os.path.join(os.path.abspath('test_models'), 'test_createPathway.sbml'))
        assert self.compareModels(inModel.getModel(), rpsbml.model)==True


    def test_createReaction(self):
        rpsbml = rpSBML.rpSBML(None, None, None, None, os.path.abspath('../cache'))
        rpsbml.genericModel('RetroPath_heterologous_pathway', 'rp_model', compXref)
        for meta in set([i for step in steps for lr in ['left', 'right'] for i in step[lr]]):
            try:
                inchi = rp_inchi[meta]
            except KeyError:
                inchi = None
            try:
                smiles = rp_smiles[meta]
            except KeyError:
                smiles = None
            rpsbml.createSpecies(meta, chemXref, None, inchi, smiles, 'MNXC3', 0, '')
        rpsbml.createPathway(path_id)
        #reactions
        step_id = 0
        for stepNum in range(len(steps)):
            rpsbml.createReaction('RP_'+str(stepNum), 'B_999999', 'B__999999', steps[stepNum], reaction_smiles[stepNum], reacXref, None, None, None)
            step_id += 1
        inModel = libsbml.readSBML(os.path.join(os.path.abspath('test_models'), 'test_createReaction.sbml'))
        assert self.compareModels(inModel.getModel(), rpsbml.model)==True


    def test_createFluxObj(self):
        rpsbml = rpSBML.rpSBML(None, None, None, None, os.path.abspath('../cache'))
        rpsbml.genericModel('RetroPath_heterologous_pathway', 'rp_model', compXref)
        for meta in set([i for step in steps for lr in ['left', 'right'] for i in step[lr]]):
            try:
                inchi = rp_inchi[meta]
            except KeyError:
                inchi = None
            try:
                smiles = rp_smiles[meta]
            except KeyError:
                smiles = None
            rpsbml.createSpecies(meta, chemXref, None, inchi, smiles, 'MNXC3', 0, '')
        rpsbml.createPathway(path_id)
        step_id = 0
        for stepNum in range(len(steps)):
            rpsbml.createReaction('RP_'+str(stepNum), 'B_999999', 'B__999999', steps[stepNum], reaction_smiles[stepNum], reacXref, None, None, None)
            step_id += 1
        rpsbml.createFluxObj('flux1', 'RP_1', 2.0, True)
        inModel = libsbml.readSBML(os.path.join(os.path.abspath('test_models'), 'test_createFluxObj.sbml'))
        assert self.compareModels(inModel.getModel(), rpsbml.model)==True

        
    def test_createMergeModels(self):
        rpsbml = rpSBML.rpSBML(None, None, None, None, os.path.abspath('../cache'))
        rpsbml.genericModel('RetroPath_heterologous_pathway', 'rp_model', compXref)
        for meta in set([i for step in steps for lr in ['left', 'right'] for i in step[lr]]):
            try:
                inchi = rp_inchi[meta]
            except KeyError:
                inchi = None
            try:
                smiles = rp_smiles[meta]
            except KeyError:
                smiles = None
            rpsbml.createSpecies(meta, chemXref, None, inchi, smiles, 'MNXC3', 0, '', None, None, None)
        rp_pathway = rpsbml.createPathway(path_id)
        #reactions
        step_id = 0
        for stepNum in range(len(steps)):
            rpsbml.createReaction('RP_'+str(stepNum), 'B_999999', 'B__999999', steps[stepNum], reaction_smiles[stepNum], reacXref, None, rpsbml.hetero_group, None)
            step_id += 1
        #other model
        document = libsbml.readSBML(os.path.join(os.path.abspath('test_models'), 'bigg_iMM904.COBRA-sbml3.xml'))
        model = document.getModel()
        rpsbml.mergeModels(model)
        inModel = libsbml.readSBML(os.path.join(os.path.abspath('test_models'), 'test_mergeModels.sbml'))
        assert self.compareModels(inModel.getModel(), model)==True
 

    def test_readSBML(self):
        rpsbml = rpSBML.rpSBML(None, None, None, None, os.path.abspath('../cache'))
        rpsbml.readSBML(os.path.join(os.path.abspath('test_models'), 'test_createReaction.sbml')) 
        inModel = libsbml.readSBML(os.path.join(os.path.abspath('test_models'), 'test_createReaction.sbml'))
        assert self.compareModels(inModel.getModel(), rpsbml.model)==True


    def test_writeSBML(self):
        rpsbml = rpSBML.rpSBML(None, None, None, None, os.path.abspath('../cache'))
        rpsbml.genericModel('RetroPath_heterologous_pathway', 'rp_model', compXref)
        rpsbml.modelName = 'test_writeSBML'
        rpsbml.writeSBML(os.path.abspath('test_models'))
        inModel = libsbml.readSBML(os.path.join(os.path.abspath('test_models'), 'test_writeSBML.sbml'))
        assert self.compareModels(inModel.getModel(), rpsbml.model)==True


    #not sure how to test this function - perhaps throw an error at it
    #def test_checklibSBML(self):

    def test_nameToSbmlId(self):
        rpsbml = rpSBML.rpSBML(None, None, None, None, os.path.abspath('../cache'))
        id_name = rpsbml._nameToSbmlId('###test_input-23####')
        assert id_name=='___test_input_23___'


    def test_genMetaID(self):
        rpsbml = rpSBML.rpSBML(None, None, None, None, os.path.abspath('../cache'))
        id_name = rpsbml._genMetaID('###test_input-23####')
        assert id_name=='_5df40e51a5d358ecfdf0372317853a79'


    def test_readAnnotation(self):
        rpsbml = rpSBML.rpSBML(None, None, None, None, os.path.abspath('../cache'))
        rpsbml.genericModel('RetroPath_heterologous_pathway', 'rp_model', compXref)
        for meta in set([i for step in steps for lr in ['left', 'right'] for i in step[lr]]):
            try:
                inchi = rp_inchi[meta]
            except KeyError:
                inchi = None
            try:
                smiles = rp_smiles[meta]
            except KeyError:
                smiles = None
            rpsbml.createSpecies(meta, chemXref, None, inchi, smiles, 'MNXC3', 0, '', None, None, None)
        annot = rpsbml.readAnnotation(rpsbml.model.getSpecies('MNXM1__64__MNXC3').getAnnotation())
        assert annot=={'bigg.metabolite': ['h', 'M_h'],
                'metanetx.chemical': ['MNXM1', 'MNXM145872', 'MNXM89553'], 
                'chebi': ['CHEBI:15378', 'CHEBI:10744', 'CHEBI:13357', 'CHEBI:5584'], 
                'hmdb': ['HMDB59597'], 
                'kegg.compound': ['C00080'], 
                'metacyc': ['PROTON'], 
                'reactome.compound': ['1132304', '113529', '1470067', '156540', '163953', '193465', '194688', '2000349', '2872447', '351626', '372511', '374900', '425969', '425978', '425999', '427899', '428040', '428548', '5668577', '70106', '74722'], 
                'seed.compound': ['39', 'cpd00067']}
        

    def test_readIBISBAAnnotation(self):
        rpsbml = rpSBML.rpSBML(None, None, None, None, os.path.abspath('../cache'))
        rpsbml.genericModel('RetroPath_heterologous_pathway', 'rp_model', compXref)
        for meta in set([i for step in steps for lr in ['left', 'right'] for i in step[lr]]):
            try:
                inchi = rp_inchi[meta]
            except KeyError:
                inchi = None
            try:
                smiles = rp_smiles[meta]
            except KeyError:
                smiles = None
            rpsbml.createSpecies(meta, chemXref, None, inchi, smiles, 'MNXC3', 0, '', None, None, None)
        rpsbml.createPathway(path_id)
        #reactions
        step_id = 0
        for stepNum in range(len(steps)):
            rpsbml.createReaction('RP_'+str(stepNum), 'B_999999', 'B__999999', steps[stepNum], reaction_smiles[stepNum], 'MNXC3', None, None)
            step_id += 1
        assert rpsbml.readIBISBAAnnotation(rpsbml.model.getSpecies('MNXM1__64__MNXC3').getAnnotation())=={'smiles': '[H+]','inchi': 'InChI=1S','inchikey': '','dG_prime_o': {'units': 'kj_per_mol', 'value': ''},'dG_prime_m': {'units': 'kj_per_mol', 'value': ''},'dG_uncert': {'units': 'kj_per_mol', 'value': ''}}


    def test_compareAnnotations(self):
        rpsbml = rpSBML.rpSBML(None, None, None, None, os.path.abspath('../cache'))
        document = libsbml.readSBML('test_models/test_createPathway.sbml')
        model = document.getModel()
        assert rpsbml.compareAnnotations(model.getSpecies('MNXM1__64__MNXC3').getAnnotation(), model.getSpecies('MNXM15__64__MNXC3').getAnnotation())==False
        assert rpsbml.compareAnnotations(model.getSpecies('MNXM1__64__MNXC3').getAnnotation(), model.getSpecies('MNXM1__64__MNXC3').getAnnotation())==True






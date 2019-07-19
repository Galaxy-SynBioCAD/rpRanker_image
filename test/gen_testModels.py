import os
import sys
import pytest
import libsbml
import pickle
import gzip

sys.path.append(os.path.join(os.path.dirname('__file__'), '..'))
import rpSBML


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


class genTestFiles():
    def __init__(self, outputPath=None):
        if outputPath==None:
            self.outputPath = 'test_models/'
            if not os.path.isdir(os.getcwd()+self.outputPath):
                os.mkdir(os.getcwd()+self.outputPath)
        else:
            self.ouputPath = outputPath
        self.compXref = pickle.load(gzip.open(os.path.join(os.path.abspath('..'), 'cache/compXref.pickle.gz'), 'rb'))
        self.chemXref = pickle.load(gzip.open(os.path.join(os.path.abspath('..'), 'cache/chemXref.pickle.gz'), 'rb'))
        self.reacXref = pickle.load(gzip.open(os.path.join(os.path.abspath('..'), 'cache/reacXref.pickle.gz'), 'rb'))

    def createModel(self):
        rpsbml = rpSBML.rpSBML('test', None, '../cache')
        rpsbml.createModel('RetroPath_heterologous_pathway', 'rp_model')
        libsbml.writeSBML(rpsbml.document, self.outputPath+'test_createModel.sbml')


    def createGenericModel(self):
        rpsbml = rpSBML.rpSBML('test', None, '../cache')
        rpsbml.genericModel('RetroPath_heterologous_pathway', 'rp_model', self.compXref)
        libsbml.writeSBML(rpsbml.document, self.outputPath+'test_genericModel.sbml')


    def createUnitDefinition(self):
        rpsbml = rpSBML.rpSBML('test', None, '../cache')
        rpsbml.createModel('RetroPath_heterologous_pathway', 'rp_model')
        unitDef = rpsbml.createUnitDefinition('mmol_per_gDW_per_hr')
        libsbml.writeSBML(rpsbml.document, self.outputPath+'test_createUnitDefinition.sbml')


    def createUnit(self):
        rpsbml = rpSBML.rpSBML('test', None, '../cache')
        rpsbml.createModel('RetroPath_heterologous_pathway', 'rp_model')
        unitDef = rpsbml.createUnitDefinition('mmol_per_gDW_per_hr')
        rpsbml.createUnit(unitDef, libsbml.UNIT_KIND_MOLE, 1, -3, 1)
        rpsbml.createUnit(unitDef, libsbml.UNIT_KIND_GRAM, 1, 0, 1)
        rpsbml.createUnit(unitDef, libsbml.UNIT_KIND_SECOND, 1, 0, 3600)
        libsbml.writeSBML(rpsbml.document, self.outputPath+'test_createUnit.sbml')


    def createParameter(self):
        rpsbml = rpSBML.rpSBML('test', None, '../cache')
        rpsbml.createModel('RetroPath_heterologous_pathway', 'rp_model')
        upInfParam = rpsbml.createParameter('B_999999', 999999.0, 'kj_per_mol')
        libsbml.writeSBML(rpsbml.document, self.outputPath+'test_createParameter.sbml')


    def createCompartent(self):
        rpsbml = rpSBML.rpSBML('test', None, '../cache')
        rpsbml.createModel('RetroPath_heterologous_pathway', 'rp_model')
        rpsbml.createCompartment(1, 'MNXC3', 'cytoplasm', self.compXref)
        libsbml.writeSBML(rpsbml.document, self.outputPath+'test_createCompartent.sbml')


    def genericModel(self):
        rpsbml = rpSBML.rpSBML('test', None, '../cache')
        rpsbml.genericModel('RetroPath_heterologous_pathway', 'rp_model', self.compXref)
        libsbml.writeSBML(rpsbml.document, self.outputPath+'test_genericModel.sbml')


    def createSpecies(self):
        rpsbml = rpSBML.rpSBML('test', None, '../cache')
        rpsbml.genericModel('RetroPath_heterologous_pathway', 'rp_model', self.compXref)
        for meta in set([i for step in steps for lr in ['left', 'right'] for i in step[lr]]):
            try:
                inchi = rp_inchi[meta]
            except KeyError:
                inchi = None
            try:
                smiles = rp_smiles[meta]
            except KeyError:
                smiles = None
            rpsbml.createSpecies(meta, self.chemXref, None, inchi, smiles, 'MNXC3', 0, '')
        libsbml.writeSBML(rpsbml.document, self.outputPath+'test_createSpecies.sbml')


    def createPathway(self):
        rpsbml = rpSBML.rpSBML('test', None, '../cache')
        rpsbml.genericModel('RetroPath_heterologous_pathway', 'rp_model', self.compXref)
        for meta in set([i for step in steps for lr in ['left', 'right'] for i in step[lr]]):
            try:
                inchi = rp_inchi[meta]
            except KeyError:
                inchi = None
            try:
                smiles = rp_smiles[meta]
            except KeyError:
                smiles = None
            rpsbml.createSpecies(meta, self.chemXref, None, inchi, smiles, 'MNXC3', 0, '')
        rpsbml.createPathway(path_id)
        #reactions
        step_id = 0
        for stepNum in range(len(steps)):
            rpsbml.createReaction('RP_'+str(stepNum), 'B_999999', 'B__999999', steps[stepNum], reaction_smiles[stepNum], self.reacXref, 'testRid', 0.0, None, rpsbml.hetero_group, None)
            step_id += 1
        libsbml.writeSBML(rpsbml.document, self.outputPath+'test_createPathway.sbml')


    def createReaction(self):
        rpsbml = rpSBML.rpSBML('test', None, '../cache')
        rpsbml.genericModel('RetroPath_heterologous_pathway', 'rp_model', self.compXref)
        for meta in set([i for step in steps for lr in ['left', 'right'] for i in step[lr]]):
            try:
                inchi = rp_inchi[meta]
            except KeyError:
                inchi = None
            try:
                smiles = rp_smiles[meta]
            except KeyError:
                smiles = None
            rpsbml.createSpecies(meta, self.chemXref, None, inchi, smiles, 'MNXC3', 0, '')
        rpsbml.createPathway(path_id)
        #reactions
        step_id = 0
        for stepNum in range(len(steps)):
            rpsbml.createReaction('RP_'+str(stepNum), 'B_999999', 'B__999999', steps[stepNum], reaction_smiles[stepNum], self.reacXref, 'testRid', 0.0, None, rpsbml.hetero_group, None)
            step_id += 1
        libsbml.writeSBML(rpsbml.document, self.outputPath+'test_createReaction.sbml')


    def createFluxObj(self):
        rpsbml = rpSBML.rpSBML('test', None, '../cache')
        rpsbml.genericModel('RetroPath_heterologous_pathway', 'rp_model', self.compXref)
        for meta in set([i for step in steps for lr in ['left', 'right'] for i in step[lr]]):
            try:
                inchi = rp_inchi[meta]
            except KeyError:
                inchi = None
            try:
                smiles = rp_smiles[meta]
            except KeyError:
                smiles = None
            rpsbml.createSpecies(meta, self.chemXref, None, inchi, smiles, 'MNXC3', 0, '')
        rpsbml.createPathway(path_id)
        step_id = 0
        for stepNum in range(len(steps)):
            rpsbml.createReaction('RP_'+str(stepNum), 'B_999999', 'B__999999', steps[stepNum], reaction_smiles[stepNum], self.reacXref, 'testRid', 0.0, None, rpsbml.hetero_group, None)
            step_id += 1
        rpsbml.createFluxObj('flux1', 'RP_1', 2.0, True)
        libsbml.writeSBML(rpsbml.document, self.outputPath+'test_createFluxObj.sbml')


    def createMergeModels(self):
        rpsbml = rpSBML.rpSBML('test', None, '../cache')
        rpsbml.genericModel('RetroPath_heterologous_pathway', 'rp_model', self.compXref)
        for meta in set([i for step in steps for lr in ['left', 'right'] for i in step[lr]]):
            try:
                inchi = rp_inchi[meta]
            except KeyError:
                inchi = None
            try:
                smiles = rp_smiles[meta]
            except KeyError:
                smiles = None
            rpsbml.createSpecies(meta, self.chemXref, None, inchi, smiles, 'MNXC3', 0, '')
        rp_pathway = rpsbml.createPathway(path_id)
        #reactions
        step_id = 0
        for stepNum in range(len(steps)):
            rpsbml.createReaction('RP_'+str(stepNum), 'B_999999', 'B__999999', steps[stepNum], reaction_smiles[stepNum], self.reacXref, 'testRid', 0.0, None, rpsbml.hetero_group, None)
            step_id += 1
        #other model
        document = libsbml.readSBML(self.outputPath+'bigg_iMM904.COBRA-sbml3.xml')
        model = document.getModel()
        rpsbml.mergeModels(model)
        libsbml.writeSBML(document, self.outputPath+'test_mergeModels.sbml')


if __name__ == "__main__":
    gtf = genTestFiles()
    print('################ createModel() ###############')
    gtf.createModel()
    print('################ createGenericModel() ###############')
    gtf.createGenericModel()
    print('################ createUnitDefinition() ###############')
    gtf.createUnitDefinition()
    print('################ createUnit() ###############')
    gtf.createUnit()
    print('################ createParameter() ###############')
    gtf.createParameter()
    print('################ createCompartent() ###############')
    gtf.createCompartent()
    print('################ createModel() ###############')
    gtf.genericModel()
    print('################ createSpecies() ###############')
    gtf.createSpecies()
    print('################ createPathway() ###############')
    gtf.createPathway()
    print('################ createReaction() ###############')
    gtf.createReaction()
    print('################ createFluxObj() ###############')
    gtf.createFluxObj()
    print('################ createMergeModels() ###############')
    gtf.createMergeModels()

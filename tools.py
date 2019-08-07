from .rpSBML import rpSBML
from io import BytesIO
import tarfile
import libsbml
import sys
import os

class tools:
    def __init__(self):
        self.pathId = 'rp_pathway'


    ## Function that takes into input a TAR collection of sbml's and opens them in memory
    #
    #
    def readrpSBMLtar(self, inputTar):
        rpsbml_paths = {}
        tar = tarfile.open(inputTar)
        for member in tar.getmembers():
            rpsbml_paths[member.name] = rpSBML(member.name,libsbml.readSBMLFromString(tar.extractfile(member).read().decode("utf-8")))
        return rpsbml_paths


    ## Function to write a collection of SBML's to TAR.XZ
    #
    #
    def writerpSBMLtar(self, rpsbml_paths, outTar):
        with tarfile.open(outTar, 'w:xz') as tf:
            for rpsbml_name in rpsbml_paths:
                data = libsbml.writeSBMLToString(rpsbml_paths[rpsbml_name].document).encode('utf-8')
                fiOut = BytesIO(data)
                info = tarfile.TarInfo(rpsbml_name)
                info.size = len(data)
                tf.addfile(tarinfo=info, fileobj=fiOut)


    ## given a reaction SMILES query Selenzyme for the Uniprot ID's and its associated score
    #
    #
    #def selenzymeREST(reaction_smile, url='http://selenzyme.synbiochem.co.uk/REST'):
    def selenzymeREST(self, reaction_smile, url):
        #columns = ['smarts', 'Seq. ID', 'Score', 'Organism Source', 'Description']
        r = requests.post( os.path.join(url, 'Query') , json={'smarts': reaction_smile} )
        res = json.loads( r.content.decode('utf-8') )
        uniprotID_score = {}
        if res['data'] is not None:
            val = json.loads( res['data'] )
        else:
            raise ValueError
        if 'Seq. ID' in val and len(val['Seq. ID'])>0:
            for ix in sorted(val['Seq. ID'], key=lambda z: int(z)):
                uniprotID_score[val['Seq. ID'][ix]] = val['Score'][ix]
        else:
            raise ValueError
        return uniprotID_score


    ## Extract the reaction SMILES from an SBML, query selenzyme and write the results back to the SBML
    #
    #
    def rpSelenzyme(self, rpsbml, url, pathId='rp_pathway'):
        groups = rpsbml.model.getPlugin('groups')
        rp_pathway = groups.getGroup(pathId)
        for member in rp_pathway.getListOfMembers():
            reaction = rpsbml.model.getReaction(member.getIdRef())
            annot = reaction.getAnnotation()
            bag_ibisba = annot.getChild('RDF').getChild('Ibisba').getChild('ibisba')
            bag_miriam = annot.getChild('RDF').getChild('Description').getChild('is').getChild('Bag')
            sel = bag_ibisba.getChild('selenzyme')
            if sel.toXMLString()=='':
                annot_string = '''
    <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
      <rdf:Ibisba rdf:about="toadd">
        <ibisba:ibisba xmlns:ibisba="http://ibisba.eu">
          <ibisba:selenzyme>
          </ibisba:selenzyme>
        </ibisba:ibisba>
      </rdf:Ibisba>
    </rdf:RDF>'''
                tmp_annot = libsbml.XMLNode.convertStringToXMLNode(annot_string)
                bag_ibisba.addChild(tmp_annot.getChild('Ibisba').getChild('ibisba').getChild('selenzyme'))
                sel = bag_ibisba.getChild('selenzyme')
            react_smiles = bag_ibisba.getChild('smiles').getChild(0).toString()
            try:
                uniprotID_score = selenzymeREST(react_smiles.replace('&gt;', '>'), url)
            except ValueError:
                continue
            for uniprot in uniprotID_score:
                ##### IBIBSA ###
                annot_string = '''
    <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
      <rdf:Ibisba rdf:about="toadd">
        <ibisba:ibisba xmlns:ibisba="http://ibisba.eu">
          <ibisba:selenzyme>
            <ibisba:'''+str(uniprot)+''' value="'''+str(uniprotID_score[uniprot])+'''" />
          </ibisba:selenzyme>
        </ibisba:ibisba>
      </rdf:Ibisba>
    </rdf:RDF>'''
                tmp_annot = libsbml.XMLNode.convertStringToXMLNode(annot_string)
                sel.addChild(tmp_annot.getChild('Ibisba').getChild('ibisba').getChild('selenzyme').getChild(uniprot))
                ''' to loop through all of them use this
                for i in range(sel.getNumChildren()):
                    print(sel.getChild(i).toXMLString())
                '''
                ##### MIRIAM ###
                # NO idea why I have to create such a large tmp annotation to add to a current one
                annot_string = '''
    <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
      <rdf:Description rdf:about="TOADD">
      <bqbiol:is>
      <rdf:Bag>
       <rdf:li rdf:resource="https://identifiers.org/uniprot/'''+str(uniprot)+'''" />
      </rdf:Bag>
      </bqbiol:is>
      </rdf:Description>
    </rdf:RDF>'''
                q = libsbml.XMLNode.convertStringToXMLNode(annot_string)
                bag_miriam.addChild(q.getChild('Description').getChild('is').getChild('Bag').getChild('li'))


    def runFBA(self, inSBMLtar, inModel):
        rpsbml_paths = readrpSBMLtar(inSBMLtar)
        for rpsbml_name in rpsbml_paths:
            #read the input sbml model
            input_rpsbml = rpRanker.rpSBML('inputMergeModel')
            input_rpsbml.readSBML(inModel)
            rpsbml_paths[rpsbml_name].mergeModels(input_rpsbml.model)
            rpfba = rpRanker.rpFBA(input_rpsbml)
            rpfba.allObj()
            ##### pass FBA results to the original model ####
            groups = rpfba.rpsbml.model.getPlugin('groups')
            rp_pathway = groups.getGroup('rp_pathway')
            for member in rp_pathway.getListOfMembers():
                reacFBA = rpfba.rpsbml.model.getReaction(member.getIdRef())
                reacIN = rpsbml_paths[rpsbml_name].model.getReaction(member.getIdRef())
                reacIN.setAnnotation(reacFBA.getAnnotation())
            #rpsbml_paths[rpsbml_name] = rpsbml_paths[rpsbml_name]
            input_rpsbml = None

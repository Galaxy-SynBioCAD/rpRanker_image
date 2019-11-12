#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to pass tarball of SBML to rpCofactors

"""


import rpSBML
import rpThermo
import libsbml
import io
import tarfile
import argparse


##
#
#
def runSingleSBML(rpthermo, member_name, rpsbml_string, path_id):
    #open one of the rp SBML files
    rpsbml = rpSBML.rpSBML(member_name, libsbml.readSBMLFromString(rpsbml_string))
    rpthermo.pathway_drG_prime_m(rpsbml, path_id)
    return libsbml.writeSBMLToString(rpsbml.document).encode('utf-8')


##
#
#
def runAllSBML(inputTar, outputTar, path_id):
    rpthermo = rpThermo.rpThermo()
    #loop through all of them and run FBA on them
    with tarfile.open(outputTar, 'w:xz') as tf:
        with tarfile.open(inputTar, 'r:xz') as in_tf:
            for member in in_tf.getmembers():
                if not member.name=='':
                    data = runSingleSBML(rpthermo,
                            member.name,
                            in_tf.extractfile(member).read().decode("utf-8"),
                            path_id)
                    fiOut = io.BytesIO(data)
                    info = tarfile.TarInfo(member.name)
                    info.size = len(data)
                    tf.addfile(tarinfo=info, fileobj=fiOut)


##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper to add cofactors to generate rpSBML collection')
    parser.add_argument('-inputTar', type=str)
    parser.add_argument('-path_id', type=str)
    parser.add_argument('-outputTar', type=str)
    params = parser.parse_args()
    runAllSBML(params.inputTar,
            params.outputTar,
            params.path_id)

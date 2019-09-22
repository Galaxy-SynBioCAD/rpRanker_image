#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to populate cofactors of a tarball of SBML

"""

import rpCofactors
import rpSBML
import tarfile
import io
import libsbml
import argparse


##
#
#
def runSingleSBML(rpcofactors, member_name, rpsbml_string, path_id, compartment_id):
    #open one of the rp SBML files
    rpsbml = rpSBML.rpSBML(member_name, libsbml.readSBMLFromString(rpsbml_string))
    #rpcofactors = rpRanker.rpCofactors()
    if rpcofactors.addCofactors(rpsbml, compartment_id, path_id):
        return libsbml.writeSBMLToString(rpsbml.document).encode('utf-8')
    else:
        return ''


##
#
#
def runAllSBML(inputTar, outputTar, path_id, compartment_id):
    #loop through all of them and run FBA on them
    rpcofactors = rpCofactors.rpCofactors()
    with tarfile.open(outputTar, 'w:xz') as tf:
        with tarfile.open(inputTar, 'r:xz') as in_tf:
            for member in in_tf.getmembers():
                if not member.name=='':
                    data = runSingleSBML(rpcofactors,
                            member.name,
                            in_tf.extractfile(member).read().decode("utf-8"),
                            path_id,
                            compartment_id)
                    if not data=='':
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
    parser.add_argument('-compartment_id', type=str)
    parser.add_argument('-outputTar', type=str)
    params = parser.parse_args()
    runAllSBML(params.inputTar,
            params.outputTar,
            params.path_id,
            params.compartment_id)

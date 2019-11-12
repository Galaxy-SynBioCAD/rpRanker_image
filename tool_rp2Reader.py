#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to query rpReader REST service

"""

import io
import rpReader
import libsbml
import tarfile
import argparse


## RetroPath2.0 reader for local packages
#
#
def rp2Reader(rp2paths_compounds, rp2_scope, rp2paths_outPaths, maxRuleIds, path_id, compartment_id, outputTar):
    rpreader = rpReader.rpReader()
    rpsbml_paths = rpreader.rp2ToSBML(rp2paths_compounds,
                        rp2_scope,
                        rp2paths_outPaths,
                        maxRuleIds),
                        path_id,
                        compartment_id)
    #pass the SBML results to a tar
    if rpsbml_paths=={}:
        raise KeyError #TODO change this to more appropriate error
    outputTar = io.BytesIO()
    with open(outputTar, 'w:xz') as tf:
        for rpsbml_name in rpsbml_paths:
            data = libsbml.writeSBMLToString(rpsbml_paths[rpsbml_name].document).encode('utf-8')
            fiOut = io.BytesIO(data)
            info = tarfile.TarInfo(name=rpsbml_name)
            info.size = len(data)
            tf.addfile(tarinfo=info, fileobj=fiOut)



##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper to parse RP2 to generate rpSBML collection')
    parser.add_argument('-rp2paths_compounds', type=str)
    parser.add_argument('-rp2_scope', type=str)
    parser.add_argument('-rp2paths_outPaths', type=str)
    parser.add_argument('-maxRuleIds', type=str)
    parser.add_argument('-path_id', type=str)
    parser.add_argument('-compartment_id', type=str)
    parser.add_argument('-outputTar', type=str)
    params = parser.parse_args()
    rp2Reader(params.rp2paths_compounds,
            params.rp2_scope,
            params.rp2paths_outPaths,
            params.maxRuleIds,
            params.path_id,
            params.compartment_id,
            params.outputTar)

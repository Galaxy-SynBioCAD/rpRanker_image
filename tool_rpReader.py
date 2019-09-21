#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to query rpReader REST service

"""
import requests
import argparse
import json


##
#
#
def rp2ReaderUpload(rp2paths_compounds, 
        rp2_scope, 
        rp2paths_outPaths, 
        maxRuleIds, 
        path_id, 
        compartment_id, 
        server, 
        outputTar):
    # Post request
    data = {'maxRuleIds': maxRuleIds, 'path_id': path_id, 'compartment_id': compartment_id}
    files = {'rp2paths_compounds': open(rp2paths_compounds, 'rb'), 
             'rp2paths_outPaths': open(rp2paths_outPaths, 'rb'),
             'rp2_scope': open(rp2_scope, 'rb'),
             'data': ('data.json', json.dumps(data))}
    r = requests.post(server+'/Query', files=files)
    r.raise_for_status()
    with open(outputTar, 'wb') as ot:
        ot.write(r.content)


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
    parser.add_argument('-server', type=str)
    parser.add_argument('-outputTar', type=str)
    params = parser.parse_args()
    rp2ReaderUpload(params.rp2paths_compounds, 
            params.rp2_scope,
            params.rp2paths_outPaths,
            params.maxRuleIds,
            params.path_id,
            params.compartment_id,
            params.server,
            params.outputTar) 

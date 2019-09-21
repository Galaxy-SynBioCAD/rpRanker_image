#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to query rpCofactors REST service

"""
import requests
import argparse
import json


##
#
#
def rpCofactors(inputTar, 
        path_id, 
        compartment_id, 
        server, 
        outputTar):
    # Post request
    data = {'path_id': path_id, 'compartment_id': compartment_id}
    files = {'inputTar': open(inputTar, 'rb'), 
             'data': ('data.json', json.dumps(data))}
    r = requests.post(server+'/Query', files=files)
    r.raise_for_status()
    with open(outputTar, 'wb') as ot:
        ot.write(r.content)


##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper to add cofactors to generate rpSBML collection')
    parser.add_argument('-inputTar', type=str)
    parser.add_argument('-path_id', type=str)
    parser.add_argument('-compartment_id', type=str)
    parser.add_argument('-server', type=str)
    parser.add_argument('-outputTar', type=str)
    params = parser.parse_args()
    rpCofactors(params.inputTar, 
            params.path_id,
            params.compartment_id,
            params.server,
            params.outputTar) 

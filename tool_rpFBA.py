#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to query rpFBA REST service

"""
import requests
import argparse
import json


##
#
#
def rpFBA(inputTar, 
        path_id,
        isMerge,
        inSBML,
        server, 
        outputTar):
    # Post request
    data = {'path_id': path_id, 'isMerge': isMerge}
    files = {'inputTar': open(inputTar, 'rb'),
             'inSBML': open(inSBML, 'rb'),
             'data': ('data.json', json.dumps(data))}
    r = requests.post(server+'/Query', files=files)
    r.raise_for_status()
    with open(outputTar, 'wb') as ot:
        ot.write(r.content)


##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper to calculate FBA to generate rpFBA collection')
    parser.add_argument('-inputTar', type=str)
    parser.add_argument('-inSBML', type=str)
    parser.add_argument('-isMerge', type=bool)
    parser.add_argument('-path_id', type=str)
    parser.add_argument('-server', type=str)
    parser.add_argument('-outputTar', type=str)
    params = parser.parse_args()
    rpFBA(params.inputTar, 
            params.path_id,
            params.isMerge,
            params.inSBML,
            params.server,
            params.outputTar) 

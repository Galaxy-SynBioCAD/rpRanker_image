#!/usr/bin/env python3
"""
Created on Mar 19

@author: Melchior du Lac
@description: Galaxy script to query rpReader REST service

"""
import requests
import argparse
from io import BytesIO
import os
import time
import json
import tarfile

#import shutil

##
#
#
def rp2ReaderUpload(rp2paths_compounds, 
        rp2_scope, 
        rp2paths_outPaths, 
        maxRuleIds, 
        pathId, 
        compartmentId, 
        server, 
        outputTar):
    # Post request
    '''
    tmpTar = BytesIO()
    with tarfile.open(fileobj=tmpTar, mode='w:xz') as tmpf:
        #rp2paths_compounds
        rp2paths_compounds_fi = open(rp2paths_compounds, mode='rb').read()
        tarinfo = tarfile.TarInfo(name='rp2paths_compounds.tsv')
        tarinfo.size = len(rp2paths_compounds_fi)
        tarinfo.mtime = time.time()
        tmpf.addfile(tarinfo, BytesIO(rp2paths_compounds_fi))
        #rp2paths_outPaths
        rp2paths_outPaths_fi = open(rp2paths_outPaths, mode='rb').read()
        tarinfo = tarfile.TarInfo(name='rp2paths_outPaths.csv')
        tarinfo.size = len(rp2paths_outPaths_fi)
        tarinfo.mtime = time.time()
        tmpf.addfile(tarinfo, BytesIO(rp2paths_outPaths_fi))
        #rp2_scope
        rp2_scope_fi = open(rp2_scope, mode='rb').read()
        tarinfo = tarfile.TarInfo(name='rp2paths_outPaths.csv')
        tarinfo.size = len(rp2_scope_fi)
        tarinfo.mtime = time.time()
        tmpf.addfile(tarinfo, BytesIO(rp2_scope_fi))
    tmpTar.seek(0)
    '''
    data = {'maxRuleIds': maxRuleIds, 'pathId': pathId, 'compartmentId': compartmentId}
    print(data)
    #files = {'file': open(bytes(tmpTar), 'rb'), 'data': ('data.json', json.dumps(data))}
    files = {'rp2paths_compounds': open(rp2paths_compounds, 'rb'), 
             'rp2paths_outPaths': open(rp2paths_outPaths, 'rb'),
             'rp2_scope': open(rp2_scope, 'rb'),
             'data': ('data.json', json.dumps(data))}
    r = requests.post(server+'/Query', files=files)
    r.raise_for_status()
    #file_like_object = BytesIO(r.content)
    #tar = tarfile.open(fileobj=file_like_object)
    #shutil.copy(tar, outputTar)
    outputTar = BytesIO(r.content)    


##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper to parse RP2 to generate rpSBML collection')
    parser.add_argument('-rp2paths_compounds', type=str)
    parser.add_argument('-rp2_scope', type=str)
    parser.add_argument('-rp2paths_outPaths', type=str)
    parser.add_argument('-maxRuleIds', type=str)
    parser.add_argument('-pathId', type=str)
    parser.add_argument('-compartmentId', type=str)
    parser.add_argument('-server', type=str)
    parser.add_argument('-outputTar', type=str)
    params = parser.parse_args()
    rp2ReaderUpload(params.rp2paths_compounds, 
            params.rp2_scope,
            params.rp2paths_outPaths,
            params.maxRuleIds,
            params.pathId,
            params.compartmentId,
            params.server,
            params.outputTar) 

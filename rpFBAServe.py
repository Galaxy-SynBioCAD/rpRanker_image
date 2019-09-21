#!/usr/bin/env python3

#from contextlib import closing
#import time
import libsbml
import argparse
import sys #exit using sys exit if any error is encountered
import os

import io
#import zipfile
import tarfile

import json
from datetime import datetime
from flask import Flask, request, jsonify, send_file, abort
from flask_restful import Resource, Api

import rpThermo
import rpSBML
import rpFBA

###################################################################################
###################################################################################
###################################################################################

import inspect
import traceback

from functools import wraps
from multiprocessing import Process, Queue


class Sentinel:
    pass


def processify(func):
    '''Decorator to run a function as a process.
    Be sure that every argument and the return value
    is *pickable*.
    The created process is joined, so the code does not
    run in parallel.
    '''

    def process_generator_func(q, *args, **kwargs):
        result = None
        error = None
        it = iter(func())
        while error is None and result != Sentinel:
            try:
                result = next(it)
                error = None
            except StopIteration:
                result = Sentinel
                error = None
            except Exception:
                ex_type, ex_value, tb = sys.exc_info()
                error = ex_type, ex_value, ''.join(traceback.format_tb(tb))
                result = None
            q.put((result, error))

    def process_func(q, *args, **kwargs):
        try:
            result = func(*args, **kwargs)
        except Exception:
            ex_type, ex_value, tb = sys.exc_info()
            error = ex_type, ex_value, ''.join(traceback.format_tb(tb))
            result = None
        else:
            error = None

        q.put((result, error))

    def wrap_func(*args, **kwargs):
        # register original function with different name
        # in sys.modules so it is pickable
        process_func.__name__ = func.__name__ + 'processify_func'
        setattr(sys.modules[__name__], process_func.__name__, process_func)

        q = Queue()
        p = Process(target=process_func, args=[q] + list(args), kwargs=kwargs)
        p.start()
        result, error = q.get()
        p.join()

        if error:
            ex_type, ex_value, tb_str = error
            message = '%s (in subprocess)\n%s' % (str(ex_value), tb_str)
            raise ex_type(message)

        return result

    def wrap_generator_func(*args, **kwargs):
        # register original function with different name
        # in sys.modules so it is pickable
        process_generator_func.__name__ = func.__name__ + 'processify_generator_func'
        setattr(sys.modules[__name__], process_generator_func.__name__, process_generator_func)

        q = Queue()
        p = Process(target=process_generator_func, args=[q] + list(args), kwargs=kwargs)
        p.start()

        result = None
        error = None
        while error is None:
            result, error = q.get()
            if result == Sentinel:
                break
            yield result
        p.join()

        if error:
            ex_type, ex_value, tb_str = error
            message = '%s (in subprocess)\n%s' % (str(ex_value), tb_str)
            raise ex_type(message)

    @wraps(func)
    def wrapper(*args, **kwargs):
        if inspect.isgeneratorfunction(func):
            return wrap_generator_func(*args, **kwargs)
        else:
            return wrap_func(*args, **kwargs)
    return wrapper


###########################################################
################## multiprocesses run #####################
###########################################################

#hack to stop the memory leak. Indeed it seems that looping through rpFBA and the rest causes a memory leak... According to: https://github.com/opencobra/cobrapy/issues/568 there is still memory leak issues with cobrapy. looping through hundreds of models and running FBA may be the culprit

@processify
#TODO: switch to merge input SBML with rpSBML model
def runSingleSBML(member_name, rpsbml_string, input_rpsbml_string, isMerge, path_id):
    #open one of the rp SBML files
    input_rpsbml = rpSBML.rpSBML('inputMergeModel', libsbml.readSBMLFromString(input_rpsbml_string))
    rpsbml = rpSBML.rpSBML(member_name, libsbml.readSBMLFromString(rpsbml_string))
    print(rpsbml)
    print(input_rpsbml)
    #read the input GEM sbml model
    #input_rpsbml = rpSBML.rpSBML('inputMergeModel', libsbml.readSBMLFromString(inSBML_string))
    #print(input_rpsbml)
    #print(input_rpsbml.model)
    rpsbml.mergeModels(input_rpsbml.model)
    rpfba = rpFBA.rpFBA(input_rpsbml)
    rpfba.allObj(path_id)
    if isMerge:
        ##### pass FBA results to the original model ####
        groups = rpfba.rpsbml.model.getPlugin('groups')
        rp_pathway = groups.getGroup(path_id)
        for member in rp_pathway.getListOfMembers():
            reacFBA = rpfba.rpsbml.model.getReaction(member.getIdRef())
            reacIN = rpsbml.model.getReaction(member.getIdRef())
            reacIN.setAnnotation(reacFBA.getAnnotation()) 
        return libsbml.writeSBMLToString(rpsbml.document).encode('utf-8')
    else:
        return libsbml.writeSBMLToString(input_rpsbml.document).encode('utf-8')


def runAllSBML(inputTar, inSBML_bytes, outTar, isMerge, path_id):
    #loop through all of them and run FBA on them
    inSBML_string = inSBML_bytes.read().decode("utf-8")
    with tarfile.open(fileobj=outTar, mode='w:xz') as tf:
        with tarfile.open(fileobj=inputTar, mode='r:xz') as in_tf:
            for member in in_tf.getmembers():
                if not member.name=='':
                    data = runSingleSBML(member.name, 
                                        in_tf.extractfile(member).read().decode("utf-8"),
                                        inSBML_string,
                                        isMerge, 
                                        path_id)
                    fiOut = io.BytesIO(data)
                    info = tarfile.TarInfo(member.name)
                    info.size = len(data)
                    tf.addfile(tarinfo=info, fileobj=fiOut)


#######################################################
############## REST ###################################
#######################################################


app = Flask(__name__)
api = Api(app)
#dataFolder = os.path.join( os.path.dirname(__file__),  'data' )


def stamp(data, status=1):
    appinfo = {'app': 'rpThermo', 'version': '1.0',
               'author': 'Melchior du Lac',
               'organization': 'BRS',
               'time': datetime.now().isoformat(),
               'status': status}
    out = appinfo.copy()
    out['data'] = data
    return out


class RestApp(Resource):
    """ REST App."""
    def post(self):
        return jsonify(stamp(None))
    def get(self):
        return jsonify(stamp(None))


class RestQuery(Resource):
    """ REST interface that generates the Design.
        Avoid returning numpy or pandas object in
        order to keep the client lighter.
    """
    def post(self):
        inputTar = request.files['inputTar']
        inSBML = request.files['inSBML']
        params = json.load(request.files['data'])
        #pass the files to the rpReader
        outputTar = io.BytesIO()
        runAllSBML(inputTar, inSBML, outputTar, bool(params['isMerge']), str(params['path_id']))
        ###### IMPORTANT ######
        outputTar.seek(0)
        #######################
        return send_file(outputTar, as_attachment=True, attachment_filename='rpThermo.tar', mimetype='application/x-tar')


api.add_resource(RestApp, '/REST')
api.add_resource(RestQuery, '/REST/Query')


if __name__== "__main__":
    #debug = os.getenv('USER') == 'mdulac'
    app.run(host="0.0.0.0", port=8994, debug=True, threaded=True)

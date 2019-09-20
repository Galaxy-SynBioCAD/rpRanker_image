#!/usr/bin/env python3

import os
import shutil
import json
import libsbml
from datetime import datetime
from flask import Flask, request, jsonify, send_file, abort
from flask_restful import Resource, Api
from io import BytesIO
import tarfile
import csv

import rpCofactors
import rpSBML

######## code taken from 


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


###############################################
###############################################
###############################################


@processify
def runSingleSBML(rpcofactors, member_name, rpsbml_string, pathId, compartment_id):
    #open one of the rp SBML files
    rpsbml = rpSBML.rpSBML(member_name, libsbml.readSBMLFromString(rpsbml_string))
    #rpcofactors = rpRanker.rpCofactors()
    if rpcofactors.addCofactors(rpsbml, compartment_id, pathId):
        return libsbml.writeSBMLToString(rpsbml.document).encode('utf-8')
    else:
        return ''


def runAllSBML(inputTar, outTar, pathId, compartment_id):
    #loop through all of them and run FBA on them
    with tarfile.open(outTar, 'w:xz') as tf:
        with tarfile.open(inputTar, 'r:xz') as in_tf:
            for member in in_tf.getmembers():
                if not member.name=='':
                    data = sunSingleSBML(rpcofactors,
                            member.name,
                            in_tf.extractfile(member).read().decode("utf-8"),
                            pathId,
                            compartment_id)
                    if not data=='':
                        fiOut = BytesIO(data)
                        info = tarfile.TarInfo(member.name)
                        info.size = len(data)
                        tf.addfile(tarinfo=info, fileobj=fiOut)

#######################################################




app = Flask(__name__)
api = Api(app)
#dataFolder = os.path.join( os.path.dirname(__file__),  'data' )


#TODO: test that it works well
#declare the rpReader globally to avoid reading the pickle at every instance
rpcofactors = rpCofactors.rpCofactors()


def stamp( data, status=1 ):
    appinfo = {'app': 'rpCofactors', 'version': '1.0', 
               'author': 'Melchior du Lac',
               'organization': 'BRS',
               'time': datetime.now().isoformat(), 
               'status': status}
    out = appinfo.copy()
    out['data'] = data
    return out


class RestApp( Resource ):
    """ REST App."""
    def post(self):
        return jsonify( stamp(None) )
    def get(self):
        return jsonify( stamp(None) )


class RestQuery(Resource):
    """ REST interface that generates the Design.
        Avoid returning numpy or pandas object in
        order to keep the client lighter.
    """
    def post(self):
        inTarFile = request.files['inTarFile'].read()
        params = json.load(request.files['params'])
        #pass the files to the rpReader
        outTar = BytesIO()
        runAllSBML(inputTar, outTar, params['pathId'], params['compartment_id'])
        return send_file(outTar, as_attachment=True)




api.add_resource(RestApp, '/REST')
api.add_resource(RestQuery, '/REST/Query')


if __name__== "__main__":
    debug = os.getenv('USER') == 'mdulac'
    app.run(host="0.0.0.0", port=8996, debug=debug, threaded=True)


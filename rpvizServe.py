'''
doeServe (c) University of Manchester 2018

doeServe is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  Pablo Carbonell, SYNBIOCHEM
@description: A REST service for OptDes 
'''
import os
import uuid
import shutil
import json
from datetime import datetime
from flask import Flask, request, jsonify,send_file
from flask_restful import Resource, Api
from rpviz.main import run

app = Flask(__name__)
api = Api(app)
dataFolder = os.path.join( os.path.dirname(__file__),  'data' )

def stamp( data, status=1 ):
    appinfo = {'app': 'rpviz', 'version': '1.0', 
               'author': 'Anaelle Badier, Pablo Carbonell',
               'organization': 'Synbiochem',
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


class RestQuery( Resource ):
    """ REST interface that generates the Design.
        Avoid returning numpy or pandas object in
        order to keep the client lighter.
    """
    def post(self):       
        file_upload = request.files['file']
        data_upload = request.files['data']
        data = json.load( data_upload )
        try:
            selenzyme_table = data['selenzyme_table']
        except:
            selenzyme_table = 'N'
        try:
            input_format = data['input_format']
        except:
            input_format = 'sbml'
        fid = str(uuid.uuid4())
        infile= os.path.join(dataFolder,fid+'.tar') 
        content = file_upload.read()
        open(infile, 'wb').write(content)
        oid = str(uuid.uuid4())
        outfolder = os.path.join( dataFolder, oid ) 
        os.mkdir( outfolder )
        print('done')
        outfile = run( infile, outfolder, selenzyme_table=selenzyme_table, typeformat=input_format, choice='5' )
        return send_file(outfile, as_attachment=True)
        #NOTE: would we not have an accumulation files with this method?... need a cleanup script
        ''' #The below will never be executed due to the above return  
        with open(outfile,'rb') as h:
            tar = h.read()
        data = {'tar': tar}
        os.remove(infile)
        shutil.rmtree( outfolder )
        return jsonify( stamp(data, 1) )
        '''

api.add_resource(RestApp, '/REST')
api.add_resource(RestQuery, '/REST/Query')

if __name__== "__main__":
    if not os.path.exists( dataFolder ):
        os.mkdir(dataFolder)
    debug = os.getenv('USER') == 'pablo'
    app.run(host="0.0.0.0", port=8998, debug=debug, threaded=True)


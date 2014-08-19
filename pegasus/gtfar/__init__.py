# Copyright 2007-2014 University Of Southern California
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

__author__ = 'Rajiv Mayani'

__VERSION__ = 0.1
import os
import errno
from flask import Flask
from flask.ext.cache import Cache
from flask.ext.restless import APIManager
from flask.ext.sqlalchemy import SQLAlchemy

#from pegasus.gtfar import views

#validExtensions = set(['txt', 'zip', 'tar', 'gz', 'bz'])

app = Flask(__name__, static_url_path='')

# Load configuration defaults
app.config.from_object('pegasus.gtfar.defaults')

#
# Database initialization
#
db = SQLAlchemy(app)
from pegasus.gtfar.models import *

db.create_all()

#
# cahe initialization
#
cache = Cache(app)


#def isValidFile(filename):
#    return '.' in filename and filename.rsplit('.', 1)[1] in validExtensions


def createRunDirectories(result):
    path = app.config['STORAGE_DIR'] + os.sep + str(result['id'])
    try:
        os.makedirs(path)
        os.makedirs(path + os.sep + 'input')
        os.makedirs(path + os.sep + 'output')
        os.makedirs(path + os.sep + 'submit')
        os.makedirs(path + os.sep + 'scratch')
        os.system("mv %s  %s" % (app.config['UPLOAD_FOLDER'] + os.sep + str(result['filename']), path + os.sep + 'input' + os.sep + str(result['filename'])))
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    return result

apiManager = APIManager(app, flask_sqlalchemy_db=db)
apiManager.create_api(Run,
                      methods = ["GET", "POST", "DELETE", "PUT"],
                      postprocessors = {
                          'POST' : [createRunDirectories]
                      })

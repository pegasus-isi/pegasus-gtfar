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

import shutil

from flask import Flask
from flask.ext.cache import Cache
from flask.ext.restless import APIManager
from flask.ext.sqlalchemy import SQLAlchemy

#
# Flask initialization
#

app = Flask(__name__, static_url_path='')
app.config.from_object('pegasus.gtfar.defaults')

#
# Database initialization
#

db = SQLAlchemy(app)
from pegasus.gtfar.models import *

db.create_all()

#
# Cache initialization
#

cache = Cache(app)

#
# Routes initialization
#

from pegasus.gtfar import views

#
# Flask Restless
#

apiManager = APIManager(app, flask_sqlalchemy_db=db)


def create_run_directories(result):
    path = os.path.join(app.config['STORAGE_DIR'], str(result['id']))

    try:
        os.makedirs(os.path.join(path, 'input'))
        os.makedirs(os.path.join(path, 'outputs'))
        os.makedirs(os.path.join(path, 'submit'))
        os.makedirs(os.path.join(path, 'scratch'))

        shutil.move(os.path.join(app.config['UPLOAD_FOLDER'], str(result['filename'])),
                    os.path.join(path, 'input', str(result['filename'])))

    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    return result


apiManager.create_api(Run,
                      methods=['GET', 'POST', 'DELETE', 'PUT', 'PATCH'],
                      postprocessors={
                          'POST': [create_run_directories]
                      })

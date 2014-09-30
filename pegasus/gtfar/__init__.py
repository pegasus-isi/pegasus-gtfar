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

__VERSION__ = '0.1-dev'

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
# Filters initialization
#

from pegasus.gtfar import filters

#
# Flask Restless
#

api_manager = APIManager(app, flask_sqlalchemy_db=db)

#
# S3 Utils Initialization
#

IS_S3_USED = (app.config['GTFAR_STAGING_SITE'] == 's3' or app.config['GTFAR_STORAGE_SITE'] == 's3')

if IS_S3_USED:
    from pegasus.gtfar.s3 import S3Utils
    s3 = S3Utils()
    app.config['GTFAR_S3_BUCKET'] = s3.get_bucket_name()

#
# Routes initialization
#

from pegasus.gtfar import views

#
# Error Handlers initialization
#

from pegasus.gtfar import errors


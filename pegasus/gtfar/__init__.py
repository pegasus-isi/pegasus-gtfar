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

apiManager = APIManager(app, flask_sqlalchemy_db=db)

#
# Routes initialization
#

from pegasus.gtfar import views

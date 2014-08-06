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


import sys
import os
from datetime import datetime
from flask import Flask
from flask.ext.cache import Cache
from flask.ext.restless import APIManager
from flask.ext.sqlalchemy import SQLAlchemy
from sqlalchemy import Column, Integer, Boolean, String, Text, DateTime
from werkzeug.utils import secure_filename


validExtensions = set(['txt', 'zip', 'tar', 'gz', 'bz'])

app = Flask(__name__, static_url_path='')
# Load configuration defaults
app.config.from_object("pegasus.gtfar.defaults")


db = SQLAlchemy(app)
cache = Cache(app)


def isValidFile(filename):
    return '.' in filename and filename.rsplit('.', 1)[1] in validExtensions

class Run(db.Model):
    __tablename__ = "runs"
    id = Column(Integer, primary_key=True)
    name = Column(Text)
    userName = Column(Text)
    filename = Column(Text)
    # Should this be text or a foreign key to a list of jobs? jobName = Column(Text, )
    status = Column(Integer)
    readLength = Column(Integer)
    mismatches = Column(Integer)
    trimUnmapped = Column(Boolean)
    mapFiltered = Column(Boolean)
    strandRule = Column(Text)
    email = Column(Text)
    gtf = Column(Text)
    genome = Column(Text)
    created = Column(DateTime, default = datetime.utcnow)

#class File(db.Model):
#    __tablename__ = "files"
#    id = Column(Integer, primary_key=True)
#    contents = Column(Text, unique=False)

db.create_all()

apiManager = APIManager(app, flask_sqlalchemy_db = db)
apiManager.create_api(Run, methods = ["GET", "POST", "DELETE", "PUT"])




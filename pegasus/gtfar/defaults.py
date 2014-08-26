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

__author__ = 'dcbriggs'

import os

# SERVER CONFIGURATION

SERVER_HOST = "127.0.0.1"
SERVER_PORT = 5000


# FLASK CONFIGURATION

DEBUG = False

# The secret key used by Flask to encrypt session keys
SECRET_KEY = os.urandom(24)

# The URI of the database
# SQLALCHEMY_DATABASE_URI = "mysql://pegasus:secret@127.0.0.1:3306/pegasus_service"
# SQLALCHEMY_DATABASE_URI = 'sqlite:///%s/.pegasus/workflow.db' % os.getenv('HOME')
SQLALCHEMY_DATABASE_URI = 'sqlite:///test.db'

UPLOAD_FOLDER = '/tmp'

# Cache Configuration
CACHE_TYPE = 'simple'

# SERVICE CONFIGURATION

# Path to the directory where the service stores all its files
STORAGE_DIR = os.path.join(os.getenv('HOME'), '.pegasus', 'gtfar')

BIN_DIR = '/data/workspace/gtfar/bin'



# CLIENT CONFIGURATION

# Service endpoint. This is only required if you install the service
# at a URL other than "http://SERVER_HOST:SERVER_PORT/".
ENDPOINT = None


# User credentials
USERNAME = ""
PASSWORD = ""

# Pegasus Configuration

PEGASUS_HOME = '/usr'

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
import sys

#
# Set OS automatically.
#


OS_TYPE = 'MACOSX' if sys.platform.startswith('darwin') else 'LINUX'


#
# GTFAR Configuration
#


GTFAR_HOME = '/lfs1/dcbriggs/repos/pegasus-gtfar'

GTFAR_BIN_DIR = os.path.join(GTFAR_HOME, 'bin')
GTFAR_DATA_DIR = os.path.join(GTFAR_HOME, 'data')

# The directory from which Pegasus workflow runs will be invoked
GTFAR_STORAGE_DIR = os.path.join(GTFAR_DATA_DIR, 'runs')

# Directory where GTFAR reference input files will be placed
GTFAR_REF_DIR = os.path.join(GTFAR_DATA_DIR, 'reference')
GTFAR_REF_GTF_DIR = os.path.join(GTFAR_REF_DIR, 'gtf')
GTFAR_REF_GENOME_DIR = os.path.join(GTFAR_REF_DIR, 'genome')

GTFAR_EXECUTION_SITE = 'condorpool'
GTFAR_STAGING_SITE = 'local'
GTFAR_STORAGE_SITE = GTFAR_STAGING_SITE

GTFAR_S3_BUCKET = 'pegasus-gtfar'

# A Temporary directory to upload files into.
UPLOAD_FOLDER = '/tmp'

# Divisor for the split size for the Dax to use
SPLIT_DIVISOR = 750000000 # 0.75 GB

#
# Pegasus Configuration
#


PEGASUS_HOME = '/usr'


#
# Server Configuration
#


SERVER_HOST = '127.0.0.1'
SERVER_PORT = 5000
ENDPOINT = None


#
# Flask Configuration
#


DEBUG = False

# Secret key to encrypt session data
SECRET_KEY = os.urandom(24)

# The URI of the database
SQLALCHEMY_DATABASE_URI = 'mysql://gtfar_user:pegasus123@127.0.0.1:3306/gtfar'

#
# Flask Cache Configuration
#

CACHE_TYPE = 'simple'

#
# User Authentication
#
USERNAME = ''
PASSWORD = ''

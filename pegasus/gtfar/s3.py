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

import os

import StringIO

import ConfigParser

from pegasus.gtfar import app

from boto.s3.connection import S3Connection


class S3Utils(object):
    MAX_DELETE = 999

    def __init__(self, s3_cfg=None):
        if s3_cfg is None:
            s3_cfg = os.path.join(os.path.expanduser('~'), '.s3cfg')

            if not os.path.isfile(s3_cfg):
                raise Exception('Unable to locate S3 configuration file')

        ini_str = open(s3_cfg, 'r').read()
        ini_fp = StringIO.StringIO(ini_str)
        config = ConfigParser.RawConfigParser()

        config.readfp(ini_fp)
        self._access_key = config.get('pegasus@amazon', 'access_key')
        self._secret_key = config.get('pegasus@amazon', 'secret_key')

        self._conn = S3Connection(self._access_key, self._secret_key)
        self._bucket = self._conn.get_bucket(app.config['GTFAR_S3_BUCKET'])

    def delete_run_dir(self, _id):
        prefix = 'data/runs/%s/' % _id
        self.delete_key(prefix)

    def delete_staging_dir(self, _id):
        prefix = 'data/runs/%s/scratch/' % _id
        self.delete_key(prefix)

    def delete_output_dir(self, _id):
        prefix = 'data/runs/%s/output/' % _id
        self.delete_key(prefix)

    def delete_key(self, key):

        # Is the key a directory?
        if key.endswith('/'):
            keys = self._bucket.list(key)
            dir_list = []
            del_list = []

            for s3_object in keys:
                name = s3_object.name

                # Directory contains a directory
                if name.endswith('/'):
                    dir_list.insert(0, name)
                else:
                    del_list.append(name)

                    if len(del_list) == S3Utils.MAX_DELETE:
                        self._bucket.delete_keys(del_list)
                        del_list = []

            if len(del_list) > 0:
                self._bucket.delete_keys(del_list)

            for i in range(0, len(dir_list), S3Utils.MAX_DELETE):
                self._bucket.delete_keys(dir_list[i:i + S3Utils.MAX_DELETE])

        else:
            self._bucket.delete_key(key)

    def get_index_files(self):
        prefix = 'data/index'
        return self._get_files(prefix)

    def get_output_files(self, _id):
        prefix = 'data/runs/%s/output' % _id
        return self._get_files(prefix)

    def _get_files(self, prefix):
        files_rs = self._bucket.list(prefix)
        files = []

        for key in files_rs:
            if not key.name.endswith('/'):
                files.append(key.name.replace(prefix, ''))

        return files
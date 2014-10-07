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

import uuid

import StringIO

import ConfigParser

from boto.exception import S3CreateError
from boto.s3.connection import S3Connection


class S3Utils(object):
    MAX_DELETE = 999
    GTFAR_S3_BUCKET_PREFIX = 'pegasus-gtfar'

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
        self._bucket = self._init_bucket(S3Utils.GTFAR_S3_BUCKET_PREFIX)

    def _init_bucket(self, bucket_prefix):
        buckets = self._conn.get_all_buckets()

        for bucket in buckets:
            if bucket.name.startswith(bucket_prefix):
                return self._conn.get_bucket(bucket.name)
        else:
            while True:
                try:
                    bucket_name = '%s-%s' % (bucket_prefix, str(uuid.uuid4()).split('-')[-1])
                    bucket = self._conn.create_bucket(bucket_name)
                    return bucket
                except S3CreateError:
                    pass

    def dir_exists(self, key):
        key = key if key.endswith('/') else key + '/'

        files = self._bucket.list(key)

        for file in files:
            return True

        return False

    def get_download_url(self, workflow_name, file_name):
        file_path = 'data/runs/%s/output/%s' % (workflow_name, file_name)
        key = self._bucket.get_key(file_path)

        if key:
            return key.generate_url(expires_in=60)

        return None

    def get_bucket_name(self):
        return self._bucket.name

    def delete_run_dir(self, _id):
        prefix = 'data/runs/%s/' % _id
        self._delete_key(prefix)

    def delete_staging_dir(self, _id):
        prefix = 'data/runs/%s/scratch/' % _id
        self._delete_key(prefix)

    def delete_output_dir(self, _id):
        prefix = 'data/runs/%s/output/' % _id
        self._delete_key(prefix)

    def _delete_key(self, key):

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
                l = dir_list[i:i + S3Utils.MAX_DELETE]
                self._bucket.delete_keys(l)

        else:
            self._bucket.delete_key(key)

    def get_index_files(self):
        prefix = 'data/index'
        return self._get_files(prefix)

    def get_output_files(self, _id):
        prefix = 'data/runs/%s/output' % _id
        files = self._get_files(prefix)
        return [(os.path.basename(name), size) for name, size in files if size > 0]

    def _get_files(self, prefix):
        files_rs = self._bucket.list(prefix)
        files = []

        for key in files_rs:
            if not key.name.endswith('/'):
                file_size = key.size
                files.append((key.name, file_size))

        return files

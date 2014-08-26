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
import errno

import shutil

from werkzeug import secure_filename

from flask import render_template, request, redirect, url_for, json

from pegasus.gtfar import app, apiManager
from pegasus.gtfar.models import Run, isValidFile


#
# Flask Restless
#

def create_run_directories(result):
    path = os.path.join(app.config['STORAGE_DIR'], str(result['id']))

    try:
        os.makedirs(os.path.join(path, 'input'))
        os.makedirs(os.path.join(path, 'config'))
        os.makedirs(os.path.join(path, 'output'))
        os.makedirs(os.path.join(path, 'submit'))
        os.makedirs(os.path.join(path, 'scratch'))

        shutil.move(os.path.join(app.config['UPLOAD_FOLDER'], str(result['filename'])),
                    os.path.join(path, 'input', str(result['filename'])))

    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    return result


def create_config(result):
    path = os.path.join(app.config['STORAGE_DIR'], str(result['id']))

    with open(os.path.join(path, 'config', 'pegasus.conf'), 'w') as conf:
        conf.write(render_template('pegasus/pegasus.conf', base_dir=path))

    with open(os.path.join(path, 'config', 'tc.txt'), 'w') as tc_txt:
        tc_txt.write(render_template('pegasus/tc.txt', bin_dir=app.config['BIN_DIR']))

    with open(os.path.join(path, 'config', 'sites.xml'), 'w') as sites_xml:
        sites_xml.write(render_template('pegasus/sites.xml', base_dir=path))

    return result


def generate_dax(result):
    path = os.path.join(app.config['STORAGE_DIR'], str(result['id']))

    from pegasus.gtfar.dax.dax import GTFAR

    gtfar = GTFAR(result['gtf'],
                  result['genome'],
                  result['id'],
                  result['filename'],
                  base_dir=path,
                  read_length=result['readLength'],
                  mismatches=result['mismatches'],
                  is_trim_unmapped=result['trimUnmapped'],
                  is_map_filtered=result['trimUnmapped'],
                  splice=True,
                  # TODO: Use value from DB
                  strand_rule='unstranded',
                  dax=os.path.join(path, '%d' % result['id']),
                  url='url',
                  email=result['email'],
                  splits=2)

    validation_results = gtfar.validate()

    if validation_results is True:
        gtfar.annotate()
        gtfar.option_filter()
        gtfar.iterative_map()
        gtfar.write_dax()


def plan_workflow(result):
    path = os.path.join(app.config['STORAGE_DIR'], str(result['id']))

    from pegasus.workflow.wrapper import PegasusWorkflow

    conf_file = os.path.join(path, 'config', 'pegasus.conf')
    dax_file = os.path.join(path, 'input', '%d' % result['id'])
    input_dir = os.path.join(path, 'input')
    submit_dir = os.path.join(path, 'submit')

    w = PegasusWorkflow(app.config['PEGASUS_HOME'], os.path.join(path, 'submit'))

    args = [
        ('--conf', conf_file),
        ('--input-dir', input_dir),
        ('--dir', submit_dir),
        ('--relative-dir', '.'),
        ('--dax', dax_file),
        ('--sites', 'condor_pool'),
        ('--staging-site', 'local'),
        ('--output-site', 'local')
    ]

    w.plan(args)


def run_workflow(result):
    submit_dir = os.path.join(app.config['STORAGE_DIR'], str(result['id']), 'submit')

    from pegasus.workflow.wrapper import PegasusWorkflow

    w = PegasusWorkflow(app.config['PEGASUS_HOME'], submit_dir)

    w.run()


apiManager.create_api(Run,
                      methods=['GET', 'POST', 'DELETE', 'PUT', 'PATCH'],
                      postprocessors={
                          'POST': [
                              create_run_directories,
                              create_config,
                              generate_dax,
                              plan_workflow,
                              run_workflow
                          ]
                      })

#
# Views
#

@app.route("/api/upload", methods=['POST'])
def upload():
    file = request.files['file']
    if file and isValidFile(file.filename):
        filename = secure_filename(file.filename)
        file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
        return redirect(url_for('index'))


@app.route("/")
def index():
    """
    Loads up the main page
    :return the template for the main page:
    """
    runs_prefix = '%(table)s%(prefix)s0.%(table)s%(prefix)s' % {'prefix': 'api', 'table': Run.__tablename__}

    api_links = {
        'runs': url_for(runs_prefix),
        'upload': url_for('upload')
    }

    return render_template('mainView.html', apiLinks=json.dumps(api_links))


@app.route("/tests")
def tests():
    """
    Loads up the testing environment
    :return: the template for the test page
    """
    return render_template("testRunner.html")

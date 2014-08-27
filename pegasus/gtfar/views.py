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

from flask import render_template, request, redirect, url_for, json, jsonify, send_from_directory

from pegasus.gtfar import app, apiManager
from pegasus.gtfar.dax.dax import GTFAR
from pegasus.gtfar.models import Run, isValidFile
from pegasus.workflow import wrapper

from pegasus.workflow.wrapper import PegasusWorkflow

#
# Application Initialization
#


@app.before_first_request
def before_first_request():
    loaded = os.path.join(app.config['GTFAR_DATA_DIR'], '.loaded')

    if not os.path.isfile(loaded):
        from pegasus.gtfar import db
        from pegasus.gtfar.models import ReplicaEntry, ReplicaAttribute

        session = db.session

        #
        # Register GTF file with JDBCRC
        #
        for path, dir_name, files in os.walk(app.config['GTFAR_REF_GTF_DIR']):

            for file_name in files:

                if file_name.startswith('.'):
                    continue

                replica = ReplicaEntry(file_name,
                                       'file://%s' % os.path.abspath(os.path.join(path, file_name)),
                                       'local',
                                       attributes=[ReplicaAttribute('gtf', 'true')])

                session.add(replica)

            session.commit()

        #
        # Register Genome files with JDBCRC
        #
        for path, dir_name, files in os.walk(app.config['GTFAR_REF_GENOME_DIR']):

            for file_name in files:

                if file_name.startswith('.'):
                    continue

                replica = ReplicaEntry(file_name,
                                       'file://%s' % os.path.abspath(os.path.join(path, file_name)),
                                       'local',
                                       attributes=[ReplicaAttribute('genome', 'true')])

                session.add(replica)

            session.commit()

        #
        # Create a marker file to indicate that
        # the reference GTF and GENOME files have been registered with the JDBCRC
        #
        open(loaded, 'a').close()


#
# Flask Restless
#


def create_run_directories(result):
    path = os.path.join(app.config['GTFAR_STORAGE_DIR'], str(result['id']))

    try:
        os.makedirs(os.path.join(path, 'input'))
        os.mkdir(os.path.join(path, 'config'))
        os.mkdir(os.path.join(path, 'output'))
        os.mkdir(os.path.join(path, 'submit'))
        os.mkdir(os.path.join(path, 'scratch'))

        shutil.move(os.path.join(app.config['UPLOAD_FOLDER'], str(result['filename'])),
                    os.path.join(path, 'input', str(result['filename'])))

    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    return result


def create_config(result):
    path = os.path.join(app.config['GTFAR_STORAGE_DIR'], str(result['id']))

    with open(os.path.join(path, 'config', 'pegasus.conf'), 'w') as conf:
        conf.write(render_template('pegasus/pegasus.conf', base_dir=path))

    with open(os.path.join(path, 'config', 'tc.txt'), 'w') as tc_txt:
        tc_txt.write(render_template('pegasus/tc.txt', bin_dir=app.config['GTFAR_BIN_DIR']))

    with open(os.path.join(path, 'config', 'sites.xml'), 'w') as sites_xml:
        sites_xml.write(render_template('pegasus/sites.xml', base_dir=path))

    return result


def generate_dax(result):
    path = os.path.join(app.config['GTFAR_STORAGE_DIR'], str(result['id']))

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
                  url='%s#/createRun' % url_for('index'),
                  email=result['email'],
                  splits=2)

    validation_results = gtfar.validate()

    if validation_results is True:
        gtfar.annotate()
        gtfar.option_filter()
        gtfar.iterative_map()
        gtfar.write_dax()


def plan_workflow(result):
    path = os.path.join(app.config['GTFAR_STORAGE_DIR'], str(result['id']))

    conf_file = os.path.join(path, 'config', 'pegasus.conf')
    input_dir = os.path.join(path, 'input')
    dax_file = os.path.join(path, '%d.dax' % result['id'])
    submit_dir = os.path.join(path, 'submit')

    workflow = PegasusWorkflow(app.config['PEGASUS_HOME'], os.path.join(path, 'submit'))

    args = [
        ('--conf', conf_file),
        ('--input-dir', input_dir),
        ('--dir', submit_dir),
        ('--relative-dir', '.'),
        ('--dax', dax_file),
        ('--sites', 'condorpool'),
        ('--staging-site', 'local'),
        ('--output-site', 'local')
    ]

    workflow.plan(args)


def run_workflow(result):
    submit_dir = os.path.join(app.config['GTFAR_STORAGE_DIR'], str(result['id']), 'submit')

    workflow = PegasusWorkflow(app.config['PEGASUS_HOME'], submit_dir)

    workflow.run()


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



@app.route("/api/upload", methods=["POST"])
def upload():
    file = request.files["file"]
    if file and isValidFile(file.filename):
        filename = secure_filename(file.filename)
        file.save(os.path.join(app.config["UPLOAD_FOLDER"], filename))
        return redirect(url_for("index"))




@app.route("/api/runs/<int:id>/status", methods=["GET"])
def getStatus(id):
    workflow = wrapper.PegasusWorkflow(app.config["PEGASUS_HOME"], app.config["GTFAR_STORAGE_DIR"] + os.sep + str(id) + os.sep + "submit")
    status = workflow.monitor(['-l'])
    # we have to change the state to be a basic data type
    return jsonify(status)

@app.route("/api/download/<path:file>")
def download(file):
    return send_from_directory(app.config["GTFAR_STORAGE_DIR"], file)

@app.route("/api/runs/<int:id>/outputs", methods=["GET"])
def getOutputFiles(id):
    files = {"objects" : []}
    for filename in os.listdir(app.config["GTFAR_STORAGE_DIR"] + os.sep + str(id) + os.sep + "output"):
        files["objects"].append({"name" : filename})
    return jsonify(files)

@app.route("/api/runs/<int:id>/logs", methods=["GET"])
def getLogs(id):
    pass

@app.route("/tests")
def tests():
    """
    Loads up the testing environment
    :return: the template for the test page
    """
    return render_template("testRunner.html")

@app.route("/")
def index():
    """
    Loads up the main page
    :return the template for the main page:
    """
    runs_prefix = '%(table)s%(prefix)s0.%(table)s%(prefix)s' % {'prefix': 'api', 'table': Run.__tablename__}

    api_links = {
        'runs': url_for(runs_prefix),
        'upload': url_for('upload'),
        'download' : '/api/download',
        'status' : '/status',
        'outputs' : '/outputs'
    }

    return render_template('mainView.html', apiLinks=json.dumps(api_links))
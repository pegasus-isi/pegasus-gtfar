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

import math

import shutil

import jinja2

from werkzeug import secure_filename

from flask import render_template, request, redirect, url_for, json, jsonify, send_from_directory

from flask.ext.restless import ProcessingException

from sqlalchemy.orm.exc import NoResultFound

from pegasus.gtfar import app, db, s3, api_manager, IS_S3_USED, __VERSION__

from pegasus.gtfar.dax.dax import GTFAR
from pegasus.gtfar.models import Run, isValidFile

from pegasus.workflow import wrapper
from pegasus.workflow.wrapper import PegasusWorkflow, StopException

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

strandRuleOptions = set(['unstranded', 'sense', 'anti-sense'])

validExtensions = set(['gz'])


class ValidationException(Exception):
    pass


def matchesAny(set, string):
    for item in set:
        if item == string:
            return 1
    return 0


def isValidFile(filename):
    return '.' in filename and filename.rsplit('.', 1)[1] in validExtensions


validation_error_messages = []


def validate_fields(data):
    global validation_error_messages
    validation_error_messages = []

    if not 'name' in data:
        validation_error_messages.append({'field': 'Name', 'message': 'You must provide a name for the run'})
    elif not data['name'].isalnum():
        validation_error_messages.append({'field': 'Name', 'message': 'Name must be an alphanumeric value'})

    if not 'filename' in data:
        validation_error_messages.append({'field': 'File', 'message': 'You must provide a file for the run'})
    elif not isValidFile(data['filename']):
        validation_error_messages.append({'field': 'File', 'message': 'File must be a GZIPPED file. (.gz)'})

    if not 'readLength' in data:
        validation_error_messages.append(
            {'field': 'Read Length', 'message': 'You must provide a read length for the run'})
    elif not type(data['readLength']) == int or not str(data['readLength']).isdigit() or not int(
            data['readLength']) >= 50 or not int(data['readLength']) <= 128:
        validation_error_messages.append(
            {'field': 'Read Length', 'message': 'Read length must be an integer between 50 and 128 (inclusive)'})
    elif 'genSplice' in data and data['genSplice'] == True:
        if not int(data['readLength']) == 75 or not int(data['readLength']) >= 100 or not int(
                data['readLength']) <= 128:
            validation_error_messages.append({'field': 'Read Length',
                                              'message': 'When the generate new splice candidate option is true the read length must be either 75 or between 100 and 128 (inclusive)'})

    if not 'mismatches' in data:
        validation_error_messages.append(
            {'field': 'Mismatches Allowed', 'message': 'You must provide the number of mismatches allowed for the run'})
    elif not type(data['mismatches']) == int or not str(data['mismatches']).isdigit() or not int(
            data['mismatches']) >= 0 or not int(data['mismatches']) <= 8:
        validation_error_messages.append({'field': 'Mismatches Allowed',
                                          'message': 'Mismatches allowed must be an integer between 0 and 8 (inclusive)'})

    if not 'strandRule' in data:
        validation_error_messages.append(
            {'field': 'Strand Rule', 'message': 'You must provide a strand rule for the run.'})
    elif not matchesAny(strandRuleOptions, data['strandRule']):
        validation_error_messages.append(
            {'field': 'Strand Rule', 'message': 'Strand Rule must be either "Unstranded", "Sense", or "Anti-Sense".'})

    if not 'genSplice' in data:
        validation_error_messages.append({'field': 'Generate New Splice Candidates',
                                          'message': 'You must specify whether or not you want the run to generate new splice candidates.'})
    elif not type(data['genSplice']) == bool:
        validation_error_messages.append(
            {'field': 'Generate New Splice Candidates', 'message': 'Generate New Splice Candidates must be a boolean.'})

    if not 'mapFiltered' in data:
        validation_error_messages.append({'field': 'Map Filtered',
                                          'message': 'You must specify whether or not you want the run to be map filtered.'})
    elif not type(data['mapFiltered']) == bool:
        validation_error_messages.append({'field': 'Map Filtered', 'message': 'Map Filtered must be a boolean.'})

    if 'email' in data:
        emails = data['email'].split(',')
        for email in emails:
            if not '@' in email:
                validation_error_messages.append(
                    {'field': 'Email', 'message': 'You must use a real, properly formatted email address'})

    if validation_error_messages:
        raise ProcessingException()


def create_run_directories(result):
    path = os.path.join(app.config['GTFAR_STORAGE_DIR'], str(result['id']))

    try:
        os.makedirs(os.path.join(path, 'input'))
        os.mkdir(os.path.join(path, 'config'))
        os.mkdir(os.path.join(path, 'output'))
        os.mkdir(os.path.join(path, 'submit'))
        os.mkdir(os.path.join(path, 'scratch'))
        shutil.move(os.path.join(app.config['UPLOAD_FOLDER'], str(result['uploadFolder']), str(result['filename'])),
                    os.path.join(path, 'input', str(result['filename']))) # There should always be one slash for the upload file name

    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    return result


def create_config(result):
    _id = str(result['id'])
    path = os.path.join(app.config['GTFAR_STORAGE_DIR'], str(result['id']))
    output_dir = os.path.join(path, 'output')

    with open(os.path.join(path, 'config', 'pegasus.conf'), 'w') as conf:
        conf.write(render_template('pegasus/pegasus.conf', base_dir=path, version=__VERSION__))

    with open(os.path.join(path, 'config', 'tc.txt'), 'w') as tc_txt:
        tc_txt.write(render_template('pegasus/tc.txt', bin_dir=app.config['GTFAR_BIN_DIR'], os=app.config['OS_TYPE']))

    with open(os.path.join(path, 'config', 'sites.xml'), 'w') as sites_xml:
        sites_xml.write(render_template('pegasus/sites.xml', _id=_id, base_dir=path, os=app.config['OS_TYPE']))

    with open(os.path.join(path, 'config', 'om.txt'), 'w') as om_txt:
        om_txt.write(
            render_template('pegasus/om.txt', _id=_id, data_dir=app.config['GTFAR_DATA_DIR'], output_dir=output_dir))

    with open(os.path.join(path, 'config', 'notifications.conf'), 'w') as notf_conf:
        notf_conf.write(render_template('pegasus/notifications.conf', base_dir=path))

    return result


def generate_dax(result):
    path = os.path.join(app.config['GTFAR_STORAGE_DIR'], str(result['id']))
    filesize = os.path.getsize(os.path.join(path, 'input', result['filename']))

    splits = max(1, int(math.floor(filesize / app.config['SPLIT_DIVISOR'])))

    gtfar = GTFAR(result['gtf'],
                  result['genome'],
                  result['id'],
                  result['filename'],
                  base_dir=path,
                  bin_dir=app.config['GTFAR_BIN_DIR'],
                  read_length=result['readLength'],
                  mismatches=result['mismatches'],
                  is_trim_unmapped=result['trimUnmapped'],
                  is_map_filtered=result['mapFiltered'],
                  clip_reads=result['genSplice'],
                  splice=True,
                  strand_rule=result['strandRule'],
                  dax=os.path.join(path, '%d' % result['id']),
                  url='%s#/createRun' % url_for('index', _external=True),
                  email=result['email'],
                  splits=splits)

    validation_results = gtfar.validate()

    if validation_results is True:
        gtfar.annotate()
        gtfar.option_filter()
        gtfar.iterative_map()
        gtfar.analyze()
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
        ('--relative-dir', 'scratch'),
        ('--relative-submit-dir', '.'),
        ('--dax', dax_file),
        ('--sites', app.config['GTFAR_EXECUTION_SITE']),
        ('--staging-site', app.config['GTFAR_STAGING_SITE']),
        ('--output-site', app.config['GTFAR_STORAGE_SITE'])
    ]

    workflow.plan(args)


def run_workflow(result):
    submit_dir = os.path.join(app.config['GTFAR_STORAGE_DIR'], str(result['id']), 'submit')

    workflow = PegasusWorkflow(app.config['PEGASUS_HOME'], submit_dir)

    workflow.run()


def remove_run_directories(instance_id, **kw):
    try:
        #
        # Delete from file-system
        #

        shutil.rmtree(os.path.join(app.config['GTFAR_STORAGE_DIR'], str(instance_id)))

        #
        # Delete from S3
        #

        if IS_S3_USED:
            s3.delete_run_dir(str(instance_id))

    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


api_manager.create_api(Run,
                       methods=['GET', 'POST', 'DELETE', 'PUT', 'PATCH'],
                       preprocessors={
                           'POST': [
                               validate_fields
                           ],
                           'DELETE': [
                               remove_run_directories,
                           ]
                       },
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


@app.route('/api/upload/<path:upload_folder>', methods=['POST'])
def upload(upload_folder):
    upload_file = request.files['file']
    print upload_file
    if upload_file and isValidFile(upload_file.filename):
        if not os.path.isdir(os.path.join(app.config['UPLOAD_FOLDER'], upload_folder)):
            os.mkdir(os.path.join(app.config['UPLOAD_FOLDER'], upload_folder))
        filename = secure_filename(upload_file.filename)
        upload_file.save(os.path.join(app.config['UPLOAD_FOLDER'], upload_folder, filename))
        return redirect(url_for('index'))
    else:
        return jsonify({"code" : 403, "message" : "Bad file type or no file presented"})


@app.route('/api/errors', methods=['GET'])
def get_errors():
    global validation_error_messages
    errors = {'errors': validation_error_messages}
    return jsonify(errors)


@app.route('/api/runs/<int:_id>/status', methods=['GET'])
def get_status(_id):
    workflow = wrapper.PegasusWorkflow(app.config['PEGASUS_HOME'],
                                       os.path.join(app.config['GTFAR_STORAGE_DIR'], str(_id), 'submit'))
    status = workflow.monitor(['-l'])
    # we have to change the state to be a basic data type
    return jsonify(status)


@app.route('/api/runs/<int:_id>/analyze', methods=['GET'])
def analyze(_id):
    workflow = wrapper.PegasusWorkflow(app.config['PEGASUS_HOME'],
                                       os.path.join(app.config['GTFAR_STORAGE_DIR'], str(_id), 'submit'))
    out = workflow.analyze()

    return out.read(), 200


@app.route('/api/download/<path:file_path>')
def download(file_path):
    return send_from_directory(app.config['GTFAR_STORAGE_DIR'], file_path)


@app.route('/api/runs/<int:_id>/outputs', methods=['GET'])
def get_output_files(_id):
    files = {'objects': []}

    if IS_S3_USED and app.config['GTFAR_STORAGE_SITE'] == 's3':
        output_files = s3.get_output_files(_id)

        for name, file_size in output_files:
            files['objects'].append({'name': name,
                                     'size': jinja2.Template('{{size|filesizeformat}}').render(size=file_size)})

    else:
        path = os.path.join(app.config['GTFAR_STORAGE_DIR'], str(_id), 'output')

        for filename in os.listdir(path):
            file_size = os.path.getsize(os.path.join(path, filename))
            files['objects'].append({'name': filename,
                                     'size': jinja2.Template('{{size|filesizeformat}}').render(size=file_size)})

    return jsonify(files)


@app.route('/api/runs/<int:_id>/logs', methods=['GET'])
def get_logs(_id):
    pass


@app.route('/api/runs/<int:_id>/stop')
def stop_run(_id):
    workflow = wrapper.PegasusWorkflow(app.config['PEGASUS_HOME'],
                                       os.path.join(app.config['GTFAR_STORAGE_DIR'], str(_id), 'submit'))
    response = {
        'status': 200,
        'reason': 'Run stopped successfully.'
    }
    try:
        session = db.session
        run = session.query(Run).filter(Run.id == _id).one()

        # pegasus-remove
        workflow.stop()

        # Update status in database
        run.status = 256
        session.add(run)
        session.commit()

    except NoResultFound:
        response = {
            'status': 400,
            'reason': 'Run (%d) not found'
        }
    except StopException:
        response = {
            'status': 500,
            'reason': 'Unable to stop the requested run.'
        }

    return jsonify(response), response['status']


@app.route('/tests')
def tests():
    """
    Loads up the testing environment
    :return: the template for the test page
    """
    return render_template('testRunner.html')


@app.route('/')
def index():
    """
    Loads up the main page
    :return the template for the main page:
    """
    runs_prefix = '%(table)s%(prefix)s0.%(table)s%(prefix)s' % {'prefix': 'api', 'table': Run.__tablename__}

    api_links = {
        'runs': url_for(runs_prefix),
        'upload': '/api/upload',
        'download': '/api/download',
        'errors': '/api/errors',
        'status': '/status',
        'outputs': '/outputs',
        'stop': '/stop'
    }

    return render_template('mainView.html', apiLinks=json.dumps(api_links))

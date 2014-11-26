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

import re

import math

import shutil

import jinja2

from werkzeug import secure_filename

from flask import render_template, request, redirect, url_for, json, jsonify, send_from_directory, make_response

from sqlalchemy.orm.exc import NoResultFound

from pegasus.gtfar import app, db, s3, api_manager, IS_S3_USED, __VERSION__

from pegasus.gtfar.models import Run, Status, isValidFile

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
        # Register Index Files available locally with JDBCRC
        #

        for path, dir_name, files in os.walk(app.config['GTFAR_DATA_DIR'] + '/index'):

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
        # Register Index Files from S3 with JDBCRC
        #

        if IS_S3_USED:
            index_files = s3.get_index_files()

            for index_file in index_files:
                file_name = os.path.basename(index_file[0])

                replica = ReplicaEntry(file_name,
                                       's3://pegasus@amazon/%s/%s' % (app.config['GTFAR_S3_BUCKET'], index_file[0]),
                                       's3',
                                       attributes=[ReplicaAttribute('index', 'true')])

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
    def __init__(self, errors=None):
        self._errors = errors

    @property
    def errors(self):
        return self._errors


def matchesAny(set, string):
    for item in set:
        if item == string:
            return 1
    return 0


def isValidFile(filename):
    return '.' in filename and filename.rsplit('.', 1)[1] in validExtensions


def workflow_name_exists(workflow_name):
    session = db.session

    try:
        session.query(Run).filter(Run.name == workflow_name).one()
    except NoResultFound:
        return False

    return True


def validate_fields(data):
    validation_error_messages = []

    if not 'name' in data:
        validation_error_messages.append({'field': 'Name', 'message': 'The name field is required and must be an alphanumeric value (underscores acceptable)'})
    elif not re.match(r'\w+$', data['name']):
        validation_error_messages.append({'field': 'Name', 'message': 'The name field is required and must be an alphanumeric value (underscores acceptable)'})
    elif workflow_name_exists(data['name']):
        validation_error_messages.append({'field': 'Workflow Name',
                                          'message': 'Workflow Name must be unique. Name %s is already in use'
                                                     % data['name']})
    elif os.path.isdir(os.path.join(app.config['GTFAR_STORAGE_DIR'], data['name'])):
        validation_error_messages.append({'field': 'Workflow Name',
                                          'message': 'Workflow Name must be unique. Name %s is already in use'
                                                     % data['name']})
    elif IS_S3_USED and s3.dir_exists('data/runs/%s' % data['name']):
        validation_error_messages.append({'field': 'Workflow Name',
                                          'message': 'Workflow Name must be unique. Name %s is already in use' %
                                                     data['name']})

    if not 'filename' in data:
        validation_error_messages.append({'field': 'File', 'message': 'You must provide a file for the run'})
    elif not isValidFile(data['filename']):
        validation_error_messages.append({'field': 'File', 'message': 'File must be a GZIPPED file. (.gz)'})

    if not 'readLength' in data:
        validation_error_messages.append(
            {'field': 'Read Length', 'message': 'The read length field is required and must be an integer between 50 and 128 (inclusive)'})
    elif not type(data['readLength']) == int or not str(data['readLength']).isdigit():
        validation_error_messages.append(
            {'field': 'Read Length', 'message': 'The read length field is required and must be an integer between 50 and 128 (inclusive)'})
    elif not int(data['readLength']) >= 50 and not int(data['readLength']) <= 128:
        validation_error_messages.append(
            {'field': 'Read Length', 'message': 'The read length field is required and must be an integer between 50 and 128 (inclusive)'})
    elif 'genSplice' in data and data['genSplice'] == True:
        if not data['readLength'] >= 75:
            validation_error_messages.append({'field': 'Read Length', 'message': 'When the Generate New Splice Candidates option is true, the read length must be >= 75'})

    if not 'mismatches' in data:
        validation_error_messages.append(
            {'field': 'Mismatches Allowed', 'message': 'The mismatches field is required and must be an integer between 0 and 8 (inclusive)'})
    elif not type(data['mismatches']) == int or not str(data['mismatches']).isdigit():
        validation_error_messages.append({'field': 'Mismatches Allowed', 'message': 'The mismatches field is required and must be an integer between 0 and 8 (inclusive)'})
    elif not int(data['mismatches']) >= 0 or not int(data['mismatches']) <= 8:
        validation_error_messages.append({'field': 'Mismatches Allowed', 'message': 'The mismatches field is required and must be an integer between 0 and 8 (inclusive)'})

    if not 'strandRule' in data:
        validation_error_messages.append(
            {'field': 'Strand Rule', 'message': 'The Strand Rule field is required and must be either "Unstranded", "Sense", or "Anti-Sense".'})
    elif not matchesAny(strandRuleOptions, data['strandRule']):
        validation_error_messages.append(
            {'field': 'Strand Rule', 'message': 'The Strand Rule field is required and must be either "Unstranded", "Sense", or "Anti-Sense".'})

    if not 'genSplice' in data:
        validation_error_messages.append({'field': 'Generate New Splice Candidates',
                                          'message': 'The Generate New Splice Candidates field is required and must be a boolean value.'})
    elif not type(data['genSplice']) == bool:
        validation_error_messages.append(
            {'field': 'Generate New Splice Candidates', 'message': 'The Generate New Splice Candidates field is required and must be a boolean value.'})

    if not 'mapFiltered' in data:
        validation_error_messages.append({'field': 'Map Filtered',
                                          'message': 'The Map Filtered field is required and must be a boolean value.'})
    elif not type(data['mapFiltered']) == bool:
        validation_error_messages.append({'field': 'Map Filtered', 'message': 'The Map Filtered field is required and must be a boolean value.'})

    if 'email' in data:
        emails = data['email'].split(',')
        for email in emails:
            if not email == "" and not '@' in email:
                validation_error_messages.append(
                    {'field': 'Email', 'message': 'Emails must be properly formatted example@domain.com'})

    if validation_error_messages:
        raise ValidationException(errors=validation_error_messages)


def create_run_directories(result):
    path = os.path.join(app.config['GTFAR_STORAGE_DIR'], str(result['name']))

    try:
        os.makedirs(os.path.join(path, 'input'))
        os.mkdir(os.path.join(path, 'config'))
        os.mkdir(os.path.join(path, 'output'))
        os.mkdir(os.path.join(path, 'submit'))
        os.mkdir(os.path.join(path, 'scratch'))
        shutil.move(os.path.join(app.config['UPLOAD_FOLDER'], str(result['uploadFolder']), str(secure_filename(result['filename']))),
                    os.path.join(path, 'input', str(secure_filename(result['filename'])))) # There should always be one slash for the upload file name

    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    return result


def create_config(result):
    _id = result['name']
    path = os.path.join(app.config['GTFAR_STORAGE_DIR'], str(result['name']))
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
    from pegasus.gtfar.dax.dax import GTFAR

    path = os.path.join(app.config['GTFAR_STORAGE_DIR'], str(result['name']))
    filesize = os.path.getsize(os.path.join(path, 'input', secure_filename(result['filename'])))

    splits = max(1, int(math.floor(filesize / app.config['SPLIT_DIVISOR'])))

    gtfar = GTFAR(result['gtf'],
                  result['genome'],
                  result['name'],
                  secure_filename(result['filename']),
                  base_dir=path,
                  bin_dir=app.config['GTFAR_BIN_DIR'],
                  read_length=result['readLength'],
                  mismatches=result['mismatches'],
                  is_trim_unmapped=result['trimUnmapped'],
                  is_map_filtered=result['mapFiltered'],
                  clip_reads=result['genSplice'],
                  splice=True,
                  strand_rule=result['strandRule'],
                  dax=os.path.join(path, '%s' % result['name']),
                  url='%s#/runs/details/%s' % (url_for('index', _external=True), result['name']),
                  email=result['email'],
                  splits=splits)

    validation_results = gtfar.validate()

    if validation_results is True:
        gtfar.annotate()
        gtfar.option_filter()
        gtfar.clip_and_parse()
        gtfar.iterative_map()
        gtfar.analyze()
        gtfar.write_dax()


def plan_workflow(result):
    path = os.path.join(app.config['GTFAR_STORAGE_DIR'], str(result['name']))

    conf_file = os.path.join(path, 'config', 'pegasus.conf')
    input_dir = os.path.join(path, 'input')
    dax_file = os.path.join(path, '%s.dax' % result['name'])
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
    submit_dir = os.path.join(app.config['GTFAR_STORAGE_DIR'], str(result['name']), 'submit')

    workflow = PegasusWorkflow(app.config['PEGASUS_HOME'], submit_dir)

    workflow.run()


def remove_run_directories(instance_id, **kw):
    try:
        #
        # Delete from file-system
        #

        shutil.rmtree(os.path.join(app.config['GTFAR_STORAGE_DIR'], instance_id), ignore_errors=True)

        #
        # Delete from S3
        #

        if IS_S3_USED:
            s3.delete_run_dir(instance_id)

    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def remove_status(instance_id, **kw):
    session = db.session
    logs = session.query(Status).filter(Status.wf_id == instance_id)

    for status in logs:
        session.delete(status)
    else:
        session.commit()


api_manager.create_api(Run,
                       methods=['GET', 'POST', 'DELETE', 'PUT', 'PATCH'],
                       preprocessors={
                           'POST': [
                               validate_fields
                           ],
                           'DELETE': [
                               remove_run_directories,
                               remove_status
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

    if upload_file and isValidFile(upload_file.filename):
        # Remove old uploaded files
        shutil.rmtree(os.path.join(app.config['UPLOAD_FOLDER'], upload_folder), ignore_errors=True)

        if not os.path.isdir(os.path.join(app.config['UPLOAD_FOLDER'], upload_folder)):
            os.mkdir(os.path.join(app.config['UPLOAD_FOLDER'], upload_folder))

        filename = secure_filename(upload_file.filename)
        upload_file.save(os.path.join(app.config['UPLOAD_FOLDER'], upload_folder, filename))

        return redirect(url_for('index'))
    else:
        return jsonify({"code": 403, "message": "Bad file type or no file presented"}), 403


@app.route('/api/runs/<string:_id>/status', methods=['GET'])
def get_status(_id):
    workflow = wrapper.PegasusWorkflow(app.config['PEGASUS_HOME'],
                                       os.path.join(app.config['GTFAR_STORAGE_DIR'], _id, 'submit'))
    status = workflow.monitor(['-l'])
    # we have to change the state to be a basic data type
    return jsonify(status)


@app.route('/api/runs/<string:_id>/analyze', methods=['GET'])
def analyze(_id):
    workflow = wrapper.PegasusWorkflow(app.config['PEGASUS_HOME'],
                                       os.path.join(app.config['GTFAR_STORAGE_DIR'], _id, 'submit'))
    out = workflow.analyze()

    if 'attach' in request.args:
        response = make_response(out.read(), 200)
        response.headers['Content-Disposition'] = 'attachment; filename="%s.err.log"' % _id

        return response

    return out.read(), 200


@app.route('/api/runs/<string:name>/input/<string:file_name>')
def download_input(name, file_name):
    path = os.path.join(app.config['GTFAR_STORAGE_DIR'], name, 'input')
    return send_from_directory(path, secure_filename(file_name))


@app.route('/api/runs/<string:name>/output/<string:file_name>')
def download_output(name, file_name):
    if IS_S3_USED:
        redirect_to = s3.get_download_url(name, file_name)
        if redirect_to:
            return redirect(redirect_to)
        else:
            return '', 404
    else:
        path = os.path.join(app.config['GTFAR_STORAGE_DIR'], name, 'output')
        return send_from_directory(path, file_name)


@app.route('/api/runs/<string:_id>/outputs', methods=['GET'])
def get_output_files(_id):
    files = {'objects': []}

    if IS_S3_USED and app.config['GTFAR_STORAGE_SITE'] == 's3':
        output_files = s3.get_output_files(_id)

        for name, file_size in output_files:
            files['objects'].append({'name': name,
                                     'size': jinja2.Template('{{size|filesizeformat}}').render(size=file_size)})

    else:
        path = os.path.join(app.config['GTFAR_STORAGE_DIR'], _id, 'output')

        for filename in os.listdir(path):
            file_size = os.path.getsize(os.path.join(path, filename))
            files['objects'].append({'name': filename,
                                     'size': jinja2.Template('{{size|filesizeformat}}').render(size=file_size)})

    return jsonify(files)


@app.route('/api/runs/<string:_id>/logs', methods=['GET'])
def get_logs(_id):
    pass


@app.route('/api/runs/<string:_id>/stop')
def stop_run(_id):
    workflow = wrapper.PegasusWorkflow(app.config['PEGASUS_HOME'],
                                       os.path.join(app.config['GTFAR_STORAGE_DIR'], _id, 'submit'))
    response = {
        'status': 200,
        'reason': 'Run stopped successfully.'
    }
    try:
        session = db.session
        run = session.query(Run).filter(Run.name == _id).one()

        # pegasus-remove
        workflow.stop()

        # Update status in database
        run.status = 256
        session.add(run)
        session.commit()

    except NoResultFound:
        response = {
            'code': 400,
            'message': 'Run not found'
        }

    return jsonify(response), response['status']


@app.route('/tests')
def tests():
    """
    Loads up the testing environment
    :return: the template for the test page
    """
    return render_template('testRunner.html')


@app.route('/gtfar/')
def index():
    """
    Loads up the main page
    :return the template for the main page:
    """
    runs_prefix = '%(table)s%(prefix)s0.%(table)s%(prefix)s' % {'prefix': 'api', 'table': Run.__tablename__}

    api_links = {
        'runs': url_for(runs_prefix),
        'upload': '/api/upload',
        'download_input': '/input',
        'download_output': '/output',
        'status': '/status',
        'outputs': '/outputs',
        'stop': '/stop',
        'sample': url_for('static', filename='sample/sample.fastq.gz'),
        'analyze': '/analyze'
    }

    return render_template('mainView.html', apiLinks=json.dumps(api_links))


@app.route('/')
def home():
    return render_template('home.html')


@app.route('/help')
def help_page():
    return render_template('help.html')


@app.route('/contact')
def contact():
    return render_template('contact.html')

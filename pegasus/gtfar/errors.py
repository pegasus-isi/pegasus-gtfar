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

from flask import jsonify

from pegasus.gtfar import app
from pegasus.gtfar.views import ValidationException

from pegasus.workflow.wrapper import PlannerException, RunnerException
from pegasus.workflow.wrapper import MonitoringException, StopException, AnalyzerException


#
# JSON Header
#

JSON_HEADERS = {'Content-Type': 'application/json; charset=utf-8'}


#
# Error handlers
#

@app.errorhandler(ValidationException)
def validation_exception(error):
    error_response = {
        'code': 'FORM_VALIDATION_FAILED',
        'errors': error.errors
    }

    return jsonify(error_response), 400, JSON_HEADERS


@app.errorhandler(PlannerException)
def planner_exception(error):
    error_response = {
        'code': 'PLANNING_FAILED',
        'message': 'Pegasus planner failed with exitcode %d' % error.exit_code
    }

    return jsonify(error_response), 500, JSON_HEADERS


@app.errorhandler(RunnerException)
def runner_exception(error):
    error_response = {
        'code': 'RUN_FAILED',
        'message': 'Pegasus runner failed with exitcode %d' % error.exit_code
    }

    return jsonify(error_response), 500, JSON_HEADERS


@app.errorhandler(MonitoringException)
def runner_exception(error):
    error_response = {
        'code': 'MONITOR_FAILED',
        'message': 'Pegasus status failed with exitcode %d' % error.exit_code
    }

    return jsonify(error_response), 500, JSON_HEADERS


@app.errorhandler(StopException)
def stop_exception(error):
    error_response = {
        'code': 'STOP_FAILED',
        'message': 'Pegasus stop failed with exitcode %d' % error.exit_code
    }

    return jsonify(error_response), 500, JSON_HEADERS


@app.errorhandler(AnalyzerException)
def analyzer_exception(error):
    error_response = {
        'code': 'ANALYZER_FAILED',
        'message': 'Pegasus analyzer failed with exitcode %d' % 1
    }

    return jsonify(error_response), 500, JSON_HEADERS


@app.errorhandler(Exception)
def fail_safe_exception(error):
    error_response = {
        'code': 'UNKNOWN_ERROR',
        'message': 'An unknown error occurred %s' % error
    }

    return jsonify(error_response), 500, JSON_HEADERS


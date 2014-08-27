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
import re
import sys
import tempfile
import subprocess

from enum import Enum


class WorkflowStates(Enum):
    __order__ = 'new create planning running success failure'

    new = -1
    create = 0
    planning = 1
    running = 2
    failing = 3
    success = 4
    failure = 5

    def __ge__(self, other):
        if self.__class__ is other.__class__:
            return self._value_ >= other._value_

        return NotImplemented

    def __gt__(self, other):
        if self.__class__ is other.__class__:
            return self._value_ > other._value_

        return NotImplemented

    def __le__(self, other):
        if self.__class__ is other.__class__:
            return self._value_ <= other._value_

        return NotImplemented

    def __lt__(self, other):
        if self.__class__ is other.__class__:
            return self._value_ < other._value_

        return NotImplemented


class Workflow(object):
    def __init__(self, pegasus_home, dag_dir=None):
        self.pegasus_home = pegasus_home
        self.dag_dir = dag_dir

    def generate_dax(self, dax, command):
        raise NotImplementedError

    def plan(self, args=None):
        raise NotImplementedError

    def run(self):
        raise NotImplementedError

    def monitor(self, args=None):
        raise NotImplementedError

    def analyze(self, args=None):
        raise NotImplementedError

    def stop(self):
        raise NotImplementedError

    @property
    def dag_dir(self):
        return self._dag_dir

    @dag_dir.setter
    def dag_dir(self, dag_dir):
        if dag_dir is not None and not os.path.isdir(dag_dir):
            raise WorkflowException

        self._dag_dir = os.path.abspath(dag_dir) if dag_dir else None

    @property
    def pegasus_home(self):
        return self._pegasus_home

    @pegasus_home.setter
    def pegasus_home(self, pegasus_home):
        if not os.path.isdir(pegasus_home):
            raise WorkflowException

        self._pegasus_home = os.path.abspath(pegasus_home)


class PegasusWorkflow(Workflow):
    def generate_dax(self, dax, command):
        raise NotImplementedError

    def plan(self, args=None):
        if not self.__valid_args(args):
            raise PlannerException

        executable = self.__get_abs_exec_path('pegasus-plan')
        args = self.__arg_to_str(args)
        ec, out, err = self.__exec_command(executable, args)

        if ec != 0:
            raise PlannerException(exit_code=ec)

    def run(self):
        executable = self.__get_abs_exec_path('pegasus-run')
        ec, out, err = self.__exec_command(executable, [self.dag_dir])

        if ec != 0:
            raise RunnerException(exit_code=ec)

    def monitor(self, args=None):
        if not self.__valid_args(args):
            raise MonitoringException

        executable = self.__get_abs_exec_path('pegasus-status')
        args = self.__arg_to_str(args)
        args.append(self.dag_dir)
        ec, out, err = self.__exec_command(executable, args)

        if ec != 0:
            raise MonitoringException(exit_code=ec)

        status = self._parse_monitoring_output(out)

        status = {
            'progress': float(status[7]),
            'state': WorkflowStates[status[8].lower()].value,
            'failed': int(status[6])
        }

        return status

    def analyze(self, args=None):
        if not self.__valid_args(args):
            raise AnalyzerException

        executable = self.__get_abs_exec_path('pegasus-analyzer')
        args = self.__arg_to_str(args)
        args.append(self.dag_dir)

        ec, out, err = self.__exec_command(executable, args)

        if ec != 0:
            raise AnalyzerException(exit_code=ec)

    def stop(self):
        executable = self.__get_abs_exec_path('pegasus-remove')
        ec, out, err = self.__exec_command(executable, [self.dag_dir])

        if ec != 0:
            raise StopException(exit_code=ec)

    def __valid_args(self, args):
        if args is not None and not hasattr(args, '__iter__'):
            return False

        return True

    def __get_abs_exec_path(self, command):
        if command and len(command) > 0:
            executable = os.path.join(self.pegasus_home, 'bin', command)
            return executable

        raise ValueError

    def __arg_to_str(self, args):
        args_list = []

        for arg in args:
            if isinstance(arg, tuple):
                args_list.append('%s' % str(arg[0]))
                for item in arg[1:]:
                    args_list.append('%s' % str(item))
            else:
                args_list.append('%s' % str(arg))

        return args_list

    def __exec_command(self, executable, args):
        command = [executable]

        if args:
            command.extend(args)

        stdout = tempfile.SpooledTemporaryFile()
        stderr = tempfile.SpooledTemporaryFile()
        exit_code = subprocess.call(command, stdout=stdout, stderr=stderr)

        if exit_code != 0:
            sys.stderr.write(' '.join(command))

        # Temporary files are opened for read and write, so at the end the file pointer should be moved to 0 for reading
        stdout.seek(0)
        stderr.seek(0)

        return exit_code, stdout, stderr

    def _parse_monitoring_output(self, out):
        header_found = False

        for line in out:

            line = line.strip()

            if header_found:
                return re.split('\s+', line)

            if line.startswith('UNRDY') or line.startswith('UNREADY'):
                header_found = True
                continue

        return False


class WorkflowException(Exception):
    def __init__(self, exit_code=None):
        self._exit_code = exit_code


class DAXGeneratorException(WorkflowException):
    pass


class MonitoringException(WorkflowException):
    pass


class PlannerException(WorkflowException):
    pass


class RunnerException(WorkflowException):
    pass


class AnalyzerException(WorkflowException):
    pass


class StopException(Exception):
    pass

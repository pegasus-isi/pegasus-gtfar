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

import os

from setuptools import setup, find_packages


# Utility function to read the README file.
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def find_package_data(dirname):
    def find_paths(dirname):
        items = []
        for fname in os.listdir(dirname):
            path = os.path.join(dirname, fname)
            if os.path.isdir(path):
                items += find_paths(path)
            elif not path.endswith('.py') and not path.endswith('.pyc'):
                items.append(path)
        return items

    items = find_paths(dirname)
    return [os.path.relpath(path, dirname) for path in items]


setup(
    name='pegasus-gtfar',
    version='0.1',
    author='Pegasus Team',
    author_email='pegasus@isi.edu',
    description='ISI & Keck School collaboration to implement a Pegasus workflow with a GUI for the GTFAR workflow.',
    long_description=read('README.md'),
    license='Apache2',
    url='https://github.com/pegasus-isi/pegasus-gtfar',
    classifiers=[
        'Topic :: Distributed Computing',
        'Topic :: Internet :: WWW/HTTP :: Application',
        'Topic :: Scientific/Engineering :: Bio-informatics',
        'License :: OSI Approved :: Apache Software License',
        'Development Status :: 3 - Alpha',
        'Environment :: Web Environment',
        'Programming Language :: Python',
        'Framework :: Flask',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Operating System :: Unix'
    ],
    namespace_packages=['pegasus'],
    packages=find_packages(),
    zip_safe=True,
    scripts=[
        'bin/gtfar-server',
        'bin/gtfar-dax-generator',
        'bin/gtfar-email',
        'bin/option-filter.sh'
    ],
    install_requires=[
        'enum34==1.0',
        'argparse==1.2.1',
        'Flask==0.10.1',
        'WTForms==2.0.1',
        'Flask-Cache==0.13.1',
        'Flask-SQLAlchemy==1.0',
        'Werkzeug==0.9.6',
        'Jinja2==2.7.3',
        'SQLAlchemy==0.9.6',
        'MySQL-python==1.2.5'
    ],
    test_suite='pegasus.tests'
)


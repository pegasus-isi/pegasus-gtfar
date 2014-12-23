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

from collections import OrderedDict


class Species(object):
    def __init__(self, name, gtf, genome, chromosomes):
        self._name = name
        self._gtf = gtf
        self._genome = genome
        self._chromosomes = chromosomes

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    @property
    def gtf(self):
        return self._gtf

    @gtf.setter
    def gtf(self, gtf):
        self._gtf = gtf

    @property
    def genome(self):
        return self._genome

    @genome.setter
    def genome(self, genome):
        self._genome = genome

    @property
    def chromosomes(self):
        return self._chromosomes

    @chromosomes.setter
    def chromosomes(self, chromosomes):
        self._chromosomes = chromosomes


#
# Species Registry
#

species = OrderedDict([
    ('human', Species('human', 'gtf', 'genome', [str(i) for i in range(1, 23)].extend(['X', 'Y', 'R', 'M']))),
    ('mouse', Species('mouse', 'gtf', 'genome', [str(i) for i in range(1, 23)].extend(['X', 'Y', 'R', 'M']))),
    ('rhesus', Species('rhesus', 'gtf', 'genome', [str(i) for i in range(1, 23)].extend(['X', 'Y', 'R', 'M'])))
])

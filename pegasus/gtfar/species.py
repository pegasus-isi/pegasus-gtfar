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
    def __init__(self, name, chromosomes):
        self._name = name
        self._chromosomes = chromosomes

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    @property
    def chromosomes(self):
        return self._chromosomes

    @chromosomes.setter
    def chromosomes(self, chromosomes):
        self._chromosomes = chromosomes


#
# Species Helper
#

def init_species_registry():
    human_chrs = [str(i) for i in range(1, 23)]
    human_chrs.extend(['X', 'Y', 'R', 'M'])

    species = OrderedDict([
        ('human', Species('human', human_chrs)),
        ('mouse', Species('mouse', human_chrs)),
        ('rhesus', Species('rhesus', human_chrs))
    ])

    return species


#
# Species Registry
#

species = init_species_registry()

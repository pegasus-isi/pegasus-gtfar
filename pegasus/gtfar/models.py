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

from datetime import datetime

from sqlalchemy import Column, Integer, Boolean, DateTime

from pegasus.gtfar import db


validExtensions = set(['gz'])


def isValidFile(filename):
    return '.' in filename and filename.rsplit('.', 1)[1] in validExtensions


class Run(db.Model):
    __tablename__ = 'runs'
    id = Column(Integer, primary_key=True)
    name = Column(db.String(64))
    userName = Column(db.String(64))
    filename = Column(db.String(255))
    uploadFolder = Column(db.String(64))
    # Should this be text or a foreign key to a list of jobs? jobName = Column(Text, )
    status = Column(Integer)
    readLength = Column(Integer)
    mismatches = Column(Integer)
    trimUnmapped = Column(Boolean)
    mapFiltered = Column(Boolean)
    genSplice = Column(Boolean)
    strandRule = Column(db.String(15))
    email = Column(db.String(255))
    gtf = Column(db.String(64))
    genome = Column(db.String(64))
    created = Column(DateTime, default=datetime.utcnow)


class Status(db.Model):
    __tablename__ = 'status'

    id = Column(Integer, primary_key=True)
    wf_id = Column(Integer)
    job_name = Column(db.String(64))
    status = Column(Integer)
    log = Column(db.String(255))


#
# Pegasus Database bases replica catalog.
#


class ReplicaEntry(db.Model):
    __tablename__ = 'rc_lfn'
    __table_args__ = (db.UniqueConstraint('lfn', 'pfn', 'site', name='sk_rc_lfn'), {'mysql_engine': 'InnoDB'})

    id = db.Column(db.Integer, primary_key=True)
    lfn = db.Column(db.String(245), nullable=False, index=True)
    pfn = db.Column(db.String(245), nullable=False)
    site = db.Column(db.String(245), nullable=True)

    attributes = db.relationship('ReplicaAttribute', backref='rc_lfn', cascade='all')

    def __init__(self, lfn, pfn, site=None, attributes=None, id=None):
        if id:
            self.id = id

        self.lfn = lfn
        self.pfn = pfn
        self.site = site
        self.attributes = []
        self._attributes = {}

        if attributes:
            self.attributes = attributes
            for attr in attributes:
                self._attributes[attr.name] = attr

    def add_attr(self, name, value):
        self.__load_attributes_map()

        if name in self._attributes:
            self._attributes[name].value = value
        else:
            attr = ReplicaAttribute(name, value)
            self.attributes.append(attr)
            self._attributes[name] = attr

    def get_attr(self, name):
        self.__load_attributes_map()

        if name in self._attributes:
            return self._attributes[name]

        return None

    def del_attr(self, name):
        self.__load_attributes_map()

        if name in self._attributes:
            for index, attribute in enumerate(self.attributes):
                if attribute.name == name:
                    break
            else:
                return

            del self.attributes[index]
            del self._attributes[name]

    def has_attr_key(self, key):
        self.__load_attributes_map()
        return True if key in self._attributes else False

    def has_attr_item(self, name, value):
        self.__load_attributes_map()
        return True if name in self._attributes and self._attributes[name].value == value else False

    def __load_attributes_map(self):
        """
        When SQLAlchemy creates the ReplicaEntry object to load results of querying the database. The constructor does
        not get called as a result, user defined attributes do not get defined.
        """
        if hasattr(self, '_attributes'):
            return
        for attr in self.attributes:
            self._attributes[attr.name] = attr

    def __repr__(self):
        return '<rc_lfn(%d, %s, %s, %s, %d)>' % (self.id, self.lfn, self.pfn, self.site, len(self.attributes))


class ReplicaAttribute(db.Model):
    __tablename__ = 'rc_attr'
    __table_args__ = (db.UniqueConstraint('id', 'name', name='sk_rc_attr'), {'mysql_engine': 'InnoDB'})

    id = db.Column(db.Integer, db.ForeignKey(ReplicaEntry.id, ondelete='CASCADE'), primary_key=True)
    name = db.Column(db.String(64), nullable=False, index=True, primary_key=True)
    value = db.Column(db.String(255), nullable=False)

    def __init__(self, name, value):
        self.name = name
        self.value = value

    def __repr__(self):
        return '<rc_attr(%s, %s)>' % (self.name, self.value)

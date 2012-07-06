"""
This module declares and imports everyting needed to
define the database entitiy classes in :mod:`.entities`.
"""

from sqlalchemy import (
        ForeignKey,
        Table, Column,
        Boolean, Integer,
        String, Text, DateTime
)
from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.declarative import declared_attr

from rnaseqlyze.core.analysis import Mixins as AnalysisMixins
from rnaseqlyze.core.srr import Mixins as RNASeqRunMixins

class _Entity(object):

    @declared_attr
    def __tablename__(cls):
        return cls.__name__.lower()

    def __setattr__(self, name, value):

        # raise an exception when setting
        # attributes that are not db columns
        if not (name.startswith('_') or hasattr(type(self), name)):
            raise Exception("'%s' is not a declared attribute" % name)

        super(_Entity, self).__setattr__(name, value)

from sqlalchemy.ext.declarative import declarative_base
Entity = declarative_base(cls=_Entity)
del declarative_base

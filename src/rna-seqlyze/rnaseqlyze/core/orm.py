"""
Map Python Objects to Database Tables

A nice tutorial showing how everything works is `here
    <http://docs.sqlalchemy.org/en/latest/orm/tutorial.html>`_.
"""

from logging import getLogger
log = getLogger(__name__)

from sqlalchemy import ForeignKey
from sqlalchemy import Table, Column
from sqlalchemy import Boolean, Integer, String, Text, DateTime
from sqlalchemy.orm.properties import RelationshipProperty
from sqlalchemy.ext.declarative import declared_attr, declarative_base

from rnaseqlyze.core import security

class Entity(object):

    @declared_attr
    def __tablename__(cls):
        return cls.__name__.lower()

    __table_args__ = {'mysql_engine': 'InnoDB'}

    def __setattr__(self, name, value):
        # raise an exception when setting
        # attributes that are not db columns
        if not name[0] == '_':
            try:
                getattr(type(self), name)
            except AttributeError:
                raise Exception("'%s' is not a declared attribute" % name)
        object.__setattr__(self, name, value)

    def __init__(self, id=None, **kwargs):

        if id: kwargs.update(id=id)

        # set1 <= set2 is the same as set1.issubset(set2)
        if not set(kwargs.keys()) <= set(self.__class__.__dict__):
            raise Exception("bad keyword arguments: %s" % list(want - have))

        for attr in kwargs:
            setattr(self, attr, kwargs[attr])

# constructor=None means:
#  "don't add a default constructor",
#  but use the standard one (__init__ declard above)
#
Entity = declarative_base(cls=Entity, constructor=None)

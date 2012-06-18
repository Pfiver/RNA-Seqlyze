"""
Map Python Objects to Database Tables
"""

# infrastructure
################

from sqlalchemy import ForeignKey
from sqlalchemy import Table, Column
from sqlalchemy import Boolean, Integer, String, Text, DateTime
from sqlalchemy.orm import relationship, backref, validates
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
                raise Exception("'%s' is not mapped to a database column" % name)
        object.__setattr__(self, name, value)

Entity = declarative_base(cls=Entity)

# entities
##########

class Analysis(Entity):
    """
    The central entity

    The whole rnaseqlyze project basically revolves around this entity.

    It represents an analysis by a researcher, as described in the srs, feature 1.
    """

    def __init__(self, **kwargs):

        if 'id' in kwargs:
            raise Exception("bad keyword argument: id; id is auto-generated")

        want = set(kwargs.keys())
        have = set(self.__class__.__dict__)
        if not want.issubset(have):
            raise Exception("bad keyword arguments: %s" % list(want - have))

        for attr in kwargs:
            setattr(self, attr, kwargs[attr])

        if not self.creation_date:
            import datetime
            self.creation_date = datetime.datetime.utcnow()

    @validates('org_accession')
    def validate_org_accession(self, name, acc):
        security.check_accession(acc)
        return acc.upper()

    id = Column(Integer, primary_key=True)

    # Genebank/Entrez id
    org_id              = Column(Integer)
    # Genebank accession
    # used to constructshared data directory name
    org_accession       = Column(String)

    inputfile_name      = Column(String)
    inputfile_type      = Column(String)
    inputfile_header    = Column(Text)

    strandspecific      = Column(Boolean)
    pairended           = Column(Boolean)
    pairendlen          = Column(Integer)

    owner_name          = Column(String, ForeignKey('user.name'))
    owner               = relationship("User", backref=backref("analyses"))

    creation_date       = Column(DateTime)
    started             = Column(Boolean)
    finished            = Column(Boolean)

    # rnaseq_study = `backref` from RNASeqStudy
    # feature_predictions = `backref` from FeaturePredictions
    # hg_tracks = `backref` from HgTrack

class User(Entity):
    """
    Constitutes a user of this service
    """
    name            = Column(String, primary_key=True)
    # analyses      = `backref` from Analysis

    def __init__(self, name):
        self.name = name

# SRA analogons

class RNASeqStudy(Entity):
    """
    Constitues an SRA "SRP" == SRA Study

    Holds the "SRP...." identifier
    """
    srp             = Column(String, primary_key=True)
    # experiments   = `backref` from RNASeqExperiment

class RNASeqExperiment(Entity):
    """
    Constitutes an SRA "SRX" == SRA Experiment
    """
    srx         = Column(String, primary_key=True)
    srp_srp     = Column(Integer, ForeignKey('rnaseqstudy.srp'))
    srp         = relationship("RNASeqStudy", backref=backref("experiments"))
    # runs      = `backref` from RNASeqRun

class RNASeqRun(Entity):
    """
    Constitutes an SRA "SRR" == SRA Run
    """
    srr         = Column(String, primary_key=True)
    srx_srx     = Column(Integer, ForeignKey('rnaseqexperiment.srx'))
    srx         = relationship("RNASeqExperiment", backref=backref("runs"))

# predictions

class FeaturePredictions(Entity):
    """
    Holds a reference to the output of a FeaturePredictor
    """
    id          = Column(Integer, primary_key=True)

    type        = Column('type', String(50))

    __mapper_args__ = {'polymorphic_on': type}

# tracks

class HgTrack(Entity):
    """
    Holds the type and filename for a UCSC Genome Browser Track
    """
    id          = Column(Integer, primary_key=True)

    type        = Column('type', String(50))

    __mapper_args__ = {'polymorphic_on': type}

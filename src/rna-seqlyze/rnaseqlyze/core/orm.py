"""
Map Python Objects to Database Tables
"""

# infrastructure
################

from sqlalchemy import ForeignKey
from sqlalchemy import Table, Column
from sqlalchemy import Boolean, Integer, String, Text, DateTime
from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.declarative import declared_attr, declarative_base

class Entity(object):
    @declared_attr
    def __tablename__(cls):
        return cls.__name__.lower()
    __table_args__ = {'mysql_engine': 'InnoDB'}

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

    id = Column(Integer, primary_key=True)
    org_accession = Column(String) # Genebank accession
                                   # maps directly to filename
    owner_name = Column(String, ForeignKey('user.name'))
    owner = relationship("User", backref=backref("analyses"))
    inputfilename = Column(String)
    creation_date = Column(DateTime)
    started = Column(Boolean)
    finished = Column(Boolean)
    strandspecific = Column(Boolean)
    pairended = Column(Boolean)
    pairendlen = Column(Integer)

    # rnaseq_study = `backref` from RNASeqStudy
    # feature_predictions = `backref` from FeaturePredictions
    # hg_tracks = `backref` from HgTrack

class User(Entity):
    """
    Constitutes a user of this service
    """
    name = Column(String, primary_key=True)
    # analyses = `backref` from Analysis

    def __init__(self, name):
        self.name = name

# SRA analogons

class RNASeqStudy(Entity):
    """
    Constitues an SRA "SRP" == SRA Study

    Holds the "SRP...." identifier
    """
    srp = Column(String, primary_key=True)
    # experiments = `backref` from RNASeqExperiment

class RNASeqExperiment(Entity):
    """
    Constitutes an SRA "SRX" == SRA Experiment
    """
    srx = Column(String, primary_key=True)
    srp_srp = Column(Integer, ForeignKey('rnaseqstudy.srp'))
    srp = relationship("RNASeqStudy", backref=backref("experiments"))
    # runs = `backref` from RNASeqRun

class RNASeqRun(Entity):
    """
    Constitutes an SRA "SRR" == SRA Run
    """
    srr = Column(String, primary_key=True)
    srx_srx = Column(Integer, ForeignKey('rnaseqexperiment.srx'))
    srx = relationship("RNASeqExperiment", backref=backref("runs"))

# predictions

class FeaturePredictions(Entity):
    """
    Holds a reference to the output of a FeaturePredictor
    """
    id = Column(Integer, primary_key=True)

    type = Column('type', String(50))
    __mapper_args__ = {'polymorphic_on': type}

# tracks

class HgTrack(Entity):
    """
    Holds the type and filename for a UCSC Genome Browser Track
    """
    id = Column(Integer, primary_key=True)

    type = Column('type', String(50))
    __mapper_args__ = {'polymorphic_on': type}

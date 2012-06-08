"""
Map Python Objects to Database Tables
"""

# infrastructure
################

from sqlalchemy import ForeignKey
from sqlalchemy import Table, Column
from sqlalchemy import Integer, String, Text
from sqlalchemy.orm import relationship, backref
from sqlalchemy.orm import sessionmaker, scoped_session
from sqlalchemy.ext.declarative import declared_attr, declarative_base

class Entity(object):
    @declared_attr
    def __tablename__(cls):
        return cls.__name__.lower()
    __table_args__ = {'mysql_engine': 'InnoDB'}

Entity = declarative_base(cls=Entity)

DBSession = scoped_session(sessionmaker())

# entities
##########

class Analysis(Entity):
    """
    The central entity

    The whole rnaseqlyze project basically revolves around this entity.

    It represents an analysis by a researcher, as described in the srs, feature 1.
    """
    id = Column(Integer, primary_key=True)
    refseq_ns = Column(String) # RefSeq accession
                               # maps directly to filename

    # rnaseq_study = `backref` from RNASeqStudy
    # feature_predictions = `backref` from FeaturePredictions
    # hg_tracks = `backref` from HgTrack

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

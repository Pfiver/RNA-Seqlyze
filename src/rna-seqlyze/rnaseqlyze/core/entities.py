"""
SQLAlchemy Database Entities

A nice tutorial showing how everything works is `here
    <http://docs.sqlalchemy.org/en/latest/orm/tutorial.html>`_.
"""

from sqlalchemy import (
        ForeignKey,
        Table, Column,
        Boolean, Integer,
        String, Text, DateTime,
)
from sqlalchemy.orm import relationship, backref

from rnaseqlyze.core.analysis import Mixins as AnalysisMixins
from rnaseqlyze.core.srr import Mixins as RNASeqRunMixins
from rnaseqlyze.core.orm import Entity

class Analysis(AnalysisMixins, Entity):

    # The order of superclasses matters!
    # AnalysisMethods.__init__ calls Entity.__init__

    """
    The central entity.
    Represents an analysis by a researcher.
    The whole rnaseqlyze project basically revolves around this entity.
    """

    id                  = Column(Integer, primary_key=True)

    org_gid             = Column(Integer)  # Organisms Genebank/Entrez gid
    org_accession       = Column(String)   # Organisms Genebank accession number

    _inputfile_name     = Column("inputfile_name", String)
    inputfile_type      = Column(String)
    _genbankfile_name   = Column("genbankfile_name", String)

    strandspecific      = Column(Boolean)
    pairended           = Column(Boolean)
    pairendminlen       = Column(Integer)
    pairendmaxlen       = Column(Integer)
    pairendcfg          = Column(Integer)

    owner               = relationship("User", backref=backref("analyses"))
    owner_name          = Column(String, ForeignKey('user.name'))

    creation_date       = Column(DateTime)
    started             = Column(Boolean)
    finished            = Column(Boolean)
    stage               = Column(String)
    error               = Column(String)

    rnaseq_run          = relationship("RNASeqRun", backref=backref("analyses"))
    rnaseq_run_srr      = Column(String, ForeignKey('rnaseqrun.srr'))

    # ft_predictions    = `backref` from FeaturePredictions
    # hg_tracks         = `backref` from HgTrack
    # galaxy_datasets   = `backref` from GalaxyDataset

class UploadSession(Entity):
    """
    Is created when somebody uploads a file
    """
    id              = Column(Integer, primary_key=True)
    analysis_id     = Column(Integer, ForeignKey(Analysis.id))
    analysis        = relationship(Analysis, uselist=False)

class User(Entity):
    """
    Constitutes a user of this service
    """
    name            = Column(String, primary_key=True)
    # analyses      = `backref` from Analysis

    def __init__(self, name):
        self.name = name

# SRA analogons

class RNASeqStudy(Entity): # stub
    """
    Constitues an SRA "SRP" == SRA Study
    """
    srp             = Column(String, primary_key=True)
    # analyses      = `backref` from Analysis
    # experiments   = `backref` from RNASeqExperiment

class RNASeqExperiment(Entity): # stub
    """
    Constitutes an SRA "SRX" == SRA Experiment
    """
    srx         = Column(String, primary_key=True)
    srp_srp     = Column(Integer, ForeignKey(RNASeqStudy.srp))
    srp         = relationship(RNASeqStudy, backref=backref("experiments"))
    # runs      = `backref` from RNASeqRun

class RNASeqRun(RNASeqRunMixins, Entity):
    """
    Constitutes an SRA "SRR" == SRA Run
    """
    srr         = Column(String, primary_key=True)
    srx_srx     = Column(Integer, ForeignKey('rnaseqexperiment.srx'))
    srx         = relationship(RNASeqExperiment, backref=backref("runs"))

class UCSCOrganism(Entity):
    """
    Holds information about the mapping of UCSC browser "db" names to
    "gene id 'title'"s, and RefSeq Accessions.
    """
    acc         = Column(String, primary_key=True)
    db          = Column(String, unique=True)
    title       = Column(String, unique=True)

class GalaxyDataset(Entity):
    """
    Holds a mapping from an analysis to a galaxy dataset id
    """
    id          = Column(String, primary_key=True)
    analysis_id = Column(Integer, ForeignKey(Analysis.id), primary_key=True)
    analysis    = relationship(Analysis, backref=backref("galaxy_datasets"))
    type        = Column(String)
    name        = Column(String)

class StageLog(Entity):
    """
    Holds the log output of one processing stage
    """
    # the primary key could be stage/analysis_id
    # but using an id automatically adds ordering
    # which comes handy, because how to order the stages otherwise ?
    id          = Column(Integer, primary_key=True)
    stage       = Column(String)
    analysis_id = Column(Integer, ForeignKey(Analysis.id))
    analysis    = relationship(Analysis, backref=backref("stage_logs"))
    text        = Column(Text)

"""
Map Python Objects to Database Tables
"""

# good read:
#   http://docs.sqlalchemy.org/en/latest/orm/tutorial.html

from sqlalchemy import ForeignKey
from sqlalchemy import Table, Column
from sqlalchemy import Boolean, Integer, String, Text, DateTime
from sqlalchemy.orm import relationship, backref, validates
from sqlalchemy.orm.properties import RelationshipProperty
from sqlalchemy.ext.declarative import declared_attr, declarative_base

from rnaseqlyze.core import security
from rnaseqlyze.core.analysis import AnalysisMixins
from rnaseqlyze.core.sra import RNASeqRunMixins

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

Entity = declarative_base(cls=Entity)

class Analysis(Entity, AnalysisMixins):
    """
    The central entity.
    Represents an analysis by a researcher.
    The whole rnaseqlyze project basically revolves around this entity.
    """

    id                  = Column(Integer, primary_key=True)

    org_gid             = Column(Integer)  # Organisms Genebank/Entrez gid
    org_accession       = Column(String)   # Organisms Genebank accession number

    inputfile_name      = Column(String)
    inputfile_type      = Column(String)

    strandspecific      = Column(Boolean)
    pairended           = Column(Boolean)
    pairendlen          = Column(Integer)

    owner               = relationship("User", backref=backref("analyses"))
    owner_name          = Column(String, ForeignKey('user.name'))

    creation_date       = Column(DateTime)
    started             = Column(Boolean)
    finished            = Column(Boolean)

    galaxy_bam_id       = Column(String)

    rnaseq_run          = relationship("RNASeqRun", backref=backref("analyses"))

    # ft_predictions    = `backref` from FeaturePredictions
    # hg_tracks         = `backref` from HgTrack

    @validates('org_accession')
    def validate_org_accession(self, attr, acc):
        security.check_valid_filename(acc)
        return acc.upper()

    @validates('inputfile_name')
    def validate_org_accession(self, attr, name):
        security.check_valid_filename(name)
        if name.find('.') < 0:
            raise Exception("Please make sure your input file has a"
                            " (meaningful) extension, like .fastq or .sra")
        return name

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

        self.create_directories()

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
    srp_srp     = Column(Integer, ForeignKey('rnaseqstudy.srp'))
    srp         = relationship("RNASeqStudy", backref=backref("experiments"))
    # runs      = `backref` from RNASeqRun

class RNASeqRun(Entity, RNASeqRunMixins):
    """
    Constitutes an SRA "SRR" == SRA Run
    """
    srr         = Column(String, primary_key=True)
    srx_srx     = Column(Integer, ForeignKey('rnaseqexperiment.srx'))
    srx         = relationship("RNASeqExperiment", backref=backref("runs"))

    @validates('srr')
    def check_srr(self, key, srr):
        import string
        assert len(srr) == 9
        assert srr[:3] == 'SRR'
        assert set(srr[3:]) < set(string.digits)
        # ... what a powerful language python is! :-)
        # http://docs.python.org/library/stdtypes.html#set
        # http://docs.python.org/library/string.html#string-constants
        return srr.upper()

    def __init__(self, srr):
        self.srr = srr
        self.create_directories()

# predictions

class FeaturePredictions(Entity): # stub
    """
    Holds a reference to the output of a FeaturePredictor
    """
    id          = Column(Integer, primary_key=True)
    type        = Column('type', String(50))

    __mapper_args__ = {'polymorphic_on': type}

# tracks

class HgTrack(Entity): # stub
    """
    Holds the type and filename for a UCSC Genome Browser Track
    """
    id          = Column(Integer, primary_key=True)
    type        = Column('type', String(50))

    __mapper_args__ = {'polymorphic_on': type}

"""
Property getters and methods for Analysis instances
"""

import logging
log = logging.getLogger(__name__)

import os
import datetime
from urllib import quote

from sqlalchemy import ForeignKey
from sqlalchemy import Table, Column
from sqlalchemy import Boolean, Integer, String, Text, DateTime
from sqlalchemy.orm import relationship, backref, validates
from sqlalchemy.orm.properties import RelationshipProperty
from sqlalchemy.ext.declarative import declared_attr, declarative_base

import rnaseqlyze
from rnaseqlyze import galaxy
from rnaseqlyze import ucscbrowser
from rnaseqlyze.core import security
from rnaseqlyze.core.orm import Entity

class AnalysisMethods(object):
    """
    Here the various analysis configurations are handled
    as transparently as possible. The properties should be
    seasy to deal with so the worker.core code doesn't get too hairy.

    .. note::

        Weather the input is an SRR identifier or a sra/fastq file is
        distinguished by checking "self.inputfile_name == None" in at least
            - the Worker
            - the analysis.pt template
            - in this file

    .. note::

        Weather the organism input is a NCBI RefSeq accession or a genbank file
        is distinguished by checking "self.genbankfile_name == None" in at least
            - the Worker
            - the analysis.pt template
            - in this file
    """

    def __init__(self, **kwargs):
        Entity.__init__(self, **kwargs)
        self.creation_date = datetime.datetime.utcnow()

        if not self.org_accession:
            acc = os.path.basename(self.inputfile_name).rsplit(".")[0]
            log.debug("setting org_accession to %s" % acc)
            self.org_accession = acc
                    

    def create_data_dir(self):
        if not os.path.isdir(self.data_dir):
            os.makedirs(self.data_dir)

    def create_gb_data_dir(self):
        if self.gb_data_dir and not os.path.isdir(self.gb_data_dir):
            os.makedirs(self.gb_data_dir)

    # org_db and hg_url (which depends upon org_db) are not set
    # as a db attribute, so old analyses where the organism was not
    # known at creation time automatically get the right url set if the
    # organism later on becomes available in the UCSC Browser

    def get_hg_url(self, org_db):
        if not self.galaxy_hgtext:
            return
        track_url = "https://" + galaxy.hostname \
                    + galaxy.dataset_display_url_template \
                        .format(dataset=self.galaxy_hgtext.id)
        hg_url = ucscbrowser.custom_track_url + \
                    ucscbrowser.custom_track_params.format(
                        org_db=org_db, track_url=quote(track_url))
        return hg_url

    def get_galaxy_id(self, name):
        for ds in self.galaxy_datasets:
            if ds.name == name:
                return ds


from os.path import join
from rnaseqlyze import (
        analyses_path,
        shared_data_path,
)

class AnalysisProperties(object):

    """
    .. note::

        - `input` means `short reads data`
        - `genbank` means `"genome" database nucleotide sequence`
    """

# data uploaded or id specified ?
# -------------------------------

    @property
    def input_uploaded(self):
        return self.inputfile_name and True

    @property
    def genbank_uploaded(self)
        return self.genbankfile_name and True

# directories
# -----------

    @property
    def data_dir(self):
        return join(analyses_path, str(self.id))

    @property
    def input_data_dir(self):
        if self.input_uploaded:
            return self.data_dir
        else:
            return self.rnaseq_run.data_dir

    @property
    def genbank_data_dir(self):
        if self.genbank_uploaded:
            return self.data_dir
        else:
            return join(shared_data_path, self.org_accession)

# short reads files
# -----------------

    @property
    def inputfile_name(self):
        if self.input_uploaded:
            return self.inputfile_name
        else:
            return self.rnaseq_run.srr + ".sra"

    @property
    def inputfile_path(self):
        return join(self.input_data_dir, self.inputfile_name)

    @property
    def inputfile_base_name(self):
        return self.inputfile_name.rsplit(".", 1)[0]

    @property
    def inputfile_fqname(self):
        return self.inputfile_base_name + ".fastq"

    @property
    def inputfile_fqpath(self):
        return self.input_data_dir + self.inputfile_fqname

    @property
    def inputfile_header(self):
        fq_file = open(self.inputfile_fqpath)
        lines = [fq_file.readline() for i in range(4)]
        log.info("Header: %s" % lines[0])
        fq_file.close()
        return "".join(lines)

# organism files
# --------------

    @property
    def genbankfile_name(self):
        if self.gb_uploaded:
            return self.genbankfile_name
        else:
            return self.org_accession + ".gb"

    @property
    def genbankfile_path(self):
        return join(self.org_data_dir, self.ganbankfile_name)

    @property
    def genbankfile_base_name(self):
        return self.genbankfile_name.rsplit(".", 1)[0]

    @property
    def genbankfile_fa_name(self):
        return self.genbankfile_base_name + ".fa"

    @property
    def genbankfile_fa_path(self):
        return join(self.genbank_data_dir, self.genbankfile_fa_name)

    @property
    def xgenbankfile_name(self):
        return self.genbankfile_base_name + ".augmented.gb"

    @property
    def xgenbankfile_path(self):
        return join(self.data_dir, self.xgenbankfile_name)

# magic galaxy_xxx attributes
# ---------------------------

    galaxy_stuff = "hgtext bam coverage hpterms".split()

    for x in galaxy_stuff: exec \
    """
        @declared_attr
        def galaxy_%s(self):
            return relationship("GalaxyDataset",
                                 uselist=False, primaryjoin="%s")
    """ % (x,
                "and_(GalaxyDataset.type == '%s',"
                     "Analysis.id == GalaxyDataset.analysis_id)" % x)

    @validates(*("galaxy_" + x for x in galaxy_stuff))
    def set_bam(self, attr, dataset):
        dataset.type=attr[7:]
        return dataset

class AnalysisValidators(object):
    @validates('org_accession')
    def validate_org_accession(self, attr, acc):
        security.check_valid_filename(acc)
        return acc.upper()

    @validates('strandspecific', 'pairended')
    def validate_boolean(self, attr, val):
        return val and True or False

    @validates('inputfile_name', 'genbankfile_name')
    def validate_x_file_name(self, attr, name):
        if '\\' in name:
            name = name.rsplit('\\', 1)[1]
        security.check_valid_filename(name)
        if name.find('.') < 0:
            raise Exception("Please make sure your input file has a"
                            " (meaningful) extension, like .fastq or .sra")
        return name

"""
Property getters and methods for Analysis instances
"""

import logging
log = logging.getLogger(__name__)

import os

import rnaseqlyze

class AnalysisMixins(object):
    """
    Here the various analysis configurations are handled
    as transparently as possible. The properties should be
    seasy to deal with so the worker.core code doesn't get too hairy.
    """
    @property
    def data_dir(self):
        return os.path.join(rnaseqlyze.analyses_path, str(self.id))

    @property
    def gb_data_dir(self):
        acc = self.org_accession
        return os.path.join(rnaseqlyze.shared_data_path, acc)

    @property
    def input_data_dir(self):
        if self.inputfile_name:
            return self.data_dir
        else:
            return self.rnaseq_run.data_dir

    @property
    def genbankfile_path(self):
        if self.genbankfile_name:
            gb_name = self.genbankfile_name
            return os.path.join(self.data_dir, gb_name)
        else:
            return os.path.join(self.gb_data_dir, self.org_accession + ".gb")

    @property
    def inputfile_path(self):
        if self.inputfile_name: # The user uploaded a file
            return os.path.join(self.data_dir, self.inputfile_name)
        else: # The user specified an SRR id
            return self.rnaseq_run.sra_path

    @property
    def inputfile_fqname(self):
        if self.inputfile_name: # The user uploaded a file
            return self.inputfile_name.rsplit('.', 1)[0] + ".fastq"
        else: # The user specified an SRR id
            return self.rnaseq_run.srr + ".fastq"

    @property
    def inputfile_fqpath(self):
        return os.path.join(self.input_data_dir, self.inputfile_fqname)

    @property
    def inputfile_header(self):
        fq_file = open(self.inputfile_fqpath)
        lines = [fq_file.readline() for i in range(4)]
        log.info("Header: %s" % lines[0])
        return "".join(lines)

    @property
    def hg_url(self):
        from rnaseqlyze import galaxy
        track_url = "http://" + galaxy.hostname + \
                        galaxy.ucsc_bam_track_template % \
                                    dict(dataset=self.galaxy_bam_id)
        from urllib import quote
        from rnaseqlyze import ucscbrowser
        hg_url = ucscbrowser.custom_track_url % \
                        dict(db="sulSol1", track_url=quote(track_url))
        return hg_url

    def create_data_dir(self):
        if not os.path.isdir(self.data_dir):
            os.makedirs(self.data_dir)

    def create_gb_data_dir(self):
        if self.gb_data_dir and not os.path.isdir(self.gb_data_dir):
            os.makedirs(self.gb_data_dir)

"""
Property getters and methods for Analysis instances
"""

import logging
log = logging.getLogger(__name__)

import os

import rnaseqlyze

class AnalysisMixins(object):
    @property
    def data_dir(self):
        return os.path.join(rnaseqlyze.analyses_path, str(self.id))

    @property
    def gb_data_dir(self):
        acc = self.org_accession
        print rnaseqlyze.shared_data_path, acc
        return os.path.join(rnaseqlyze.shared_data_path, acc)

    @property
    def inputfile_path(self):
        return os.path.join(self.data_dir, self.inputfile_name)

    @property
    def inputfile_fqname(self):
        return self.inputfile_name.rsplit('.', 1)[0] + ".fastq"

    @property
    def inputfile_fqpath(self):
        return os.path.join(self.data_dir, self.inputfile_fqname)

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

    def create_directories(self):
        if not os.path.isdir(self.data_dir):
            os.makedirs(self.data_dir)
        if not os.path.isdir(self.gb_data_dir):
            os.makedirs(self.gb_data_dir)

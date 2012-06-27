"""
Sequence Run Archive interaction
"""

import logging
log = logging.getLogger(__name__)

import os
from os import path

import rnaseqlyze

srr_url_template = "http://ftp-private.ncbi.nlm.nih.gov" \
        "/sra/sra-instant/reads/ByRun/sra/SRR/{srr:.6}/{srr}/{srr}.sra"
# e.g.  "/sra/sra-instant/reads/ByRun/sra/SRR/SRR000/SRR000001/SRR000001.sra"

class RNASeqRunMixins(object):

    @property
    def data_dir(self):
        return path.join(rnaseqlyze.shared_data_path, self.srr)

    @property
    def sra_path(self):
        return path.join(self.data_dir, self.srr) + ".sra"

    def create_directories(self):
        if not os.path.isdir(self.data_dir):
            os.makedirs(self.data_dir)

    def download(self):
        from urllib import urlretrieve
        log.debug("fetching " + self.srr)
        urlretrieve(srr_url_template.format(srr=self.srr), self.sra_path)
        log.debug("done")

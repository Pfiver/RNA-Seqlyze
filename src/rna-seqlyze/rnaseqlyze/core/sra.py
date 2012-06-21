"""
Sequence Run Archive interaction
"""

import logging
log = logging.getLogger(__name__)

import os
from os import path

import rnaseqlyze

srr_url_template = "http://ftp-private.ncbi.nlm.nih.gov" \
        "/sra/sra-instant/reads/ByRun/sra/SRR/{srr:.3}/{srr:.6}/{srr}.sra"
# e.g.  "/sra/sra-instant/reads/ByRun/sra/SRR/SRR000/SRR000001/SRR000001.sra"

class RNASeqRunMixins(object):

    @property
    def srr_data_dir(self):
        if not self.srr:
            return
        return path.join(rnaseqlyze.shared_data_path, self.srr)

    def create_directories(self):
        if not os.path.isdir(self.srr_data_dir):
            os.makedirs(self.srr_data_dir)

    def download(self):
        from urllib import urlretrievea
        log.debug("fetching " + self.srr)
        urlretrieve(url, self.srr_data_dir)
        log.debug("done")

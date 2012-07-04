"""
Sequence Run Archive interaction
"""

import logging
log = logging.getLogger(__name__)

import os
from os import path
from urllib2 import urlopen
from shutil import copyfileobj

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
        return path.join(self.data_dir, self.sra_name)

    @property
    def sra_name(self):
        return self.srr + ".sra"

    def create_directories(self):
        if not os.path.isdir(self.data_dir):
            os.makedirs(self.data_dir)

    def download(self):
        try:
            log.debug("fetching " + self.srr)
            srr_url = srr_url_template.format(srr=self.srr)
            remote = urlopen(srr_url, timeout=60)
            local = open(self.sra_path, "w")
            copyfileobj(remote, local)
        except Exception, e:
            log.error("Error fetching SRR: %r" % e)
            os.unlink(self.sra_path)
            raise
        finally:
            # note:
            #  in case of an error, unlinking wil precede closing
            #   -- no problem on unix
            if local:
                local.close()
        log.debug("done")

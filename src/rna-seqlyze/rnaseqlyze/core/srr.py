"""
Sequence Run Archive interaction
"""

import logging
log = logging.getLogger(__name__)

import os
from os import path
from urllib2 import urlopen
from shutil import copyfileobj

from sqlalchemy.orm import validates

import rnaseqlyze

url_template = "http://ftp-private.ncbi.nlm.nih.gov" \
        "/sra/sra-instant/reads/ByRun/sra/SRR/{srr:.6}/{srr}/{srr}.sra"
# e.g.  "/sra/sra-instant/reads/ByRun/sra/SRR/SRR000/SRR000001/SRR000001.sra"

class Methods(object):
    def download(self):
        try:
            log.info("fetching " + self.srr)
            srr_url = url_template.format(srr=self.srr)
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

class Properties(object):
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

class Validators(object):
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

class Mixins(Methods, Properties, Validators):
    pass

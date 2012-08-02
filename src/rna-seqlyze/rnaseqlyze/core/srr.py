"""
Sequence Run Archive interaction
"""

import logging
log = logging.getLogger(__name__)

import os
from os import path
from time import time
from urlparse import urlparse
from httplib import HTTPConnection
from socket import timeout

from sqlalchemy.orm import validates

import rnaseqlyze

url_template = "http://ftp-private.ncbi.nlm.nih.gov" \
        "/sra/sra-instant/reads/ByRun/sra/{srr:.3}/{srr:.6}/{srr}/{srr}.sra"
# e.g.  "/sra/sra-instant/reads/ByRun/sra/SRR/SRR000/SRR000001/SRR000001.sra"

class Methods(object):

    def download(self):
        srr_url = url_template.format(srr=self.srr)
        log.info("fetching %s" % srr_url)
        url_parts = urlparse(srr_url)

        resp = None
        max_tries = 3

        for tries_left in reversed(range(max_tries)):
            try:
                conn = HTTPConnection(url_parts.netloc, strict=True, timeout=60)

                msg = "connecting to server"
                log.debug(msg)
                conn.connect()

                msg = "sending GET request"
                log.debug(msg)
                conn.request("GET", url_parts.path)

                msg = "waiting for response from server"
                log.debug(msg)
                resp = conn.getresponse()

                break
            except timeout:
                log.warn("timeout " + msg)
                if tries_left > 0:
                    log.info("trying %d more time%s" %
                                    (tries_left, tries_left > 1 and 's' or ''))
                else:
                    log.error("tried %d times - giving up", max_tries)
                    raise

        # Impossible really, but...
        assert resp != None

        if resp.status != 200:
            problem = "Bad response from server: %d - %s" % (
                                         resp.status, resp.msg.status)
            log.error(problem)
            raise Exception(problem)

        local = open(self.sra_path, "w")
        total = resp.length
        read = 0
        then = time()
        log.info("transfering %d kb data..." % (total / 1024
                                                    if total else -1))
        try:
            while True:
                buf = resp.read(16*1024)
                read += len(buf)
                now = time()
                if then < now - 15:
                    then = now
                    log.info("%d kb left" % (((total or 0) - read) / 1024))
                if not buf:
                    break
                local.write(buf)
        except Exception, e:
            log.error("Error downloading SRR: %r" % e)
            os.unlink(self.sra_path)
            raise
        finally:
            # note:
            #  in case of an error, unlinking wil precede closing
            #   -- no problem on unix
            local.close()

        log.debug("Success!")

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

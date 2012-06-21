import logging
log = logging.getLogger(__name__)

import os
import urllib2

import rnaseqlyze
from rnaseqlyze.core import security
from rnaseqlyze.core.orm import Analysis, User

def create_analysis(session, inputfile, attributes):
    # owner handling
    if 'owner' not in attributes:
        owner = session.query(User).get("anonymous")
        if not owner:
            owner = User("anonymous")
            session.add(owner)
        attributes['owner'] = owner
    # srr handling
    if 'sra_run' in attributes:
        try:
            attributes['sra_run'] = RNASeqRun(attributes['sra_run'])
            attributes['sra_run'].create_directories()
        except:
            # The RNASeqRun constructor checks if the SRAnnnnnn argument
            # and raises an exception unless it passes the checks
            # e.g. if field was left blank/at default value
            pass # TODO: decide/document what to do
    # create db object
    analysis = Analysis(**attributes)
    session.add(analysis)
    session.flush() # sets analysis.id
    analysis.create_directories()
    # transfer inputfile
    if self.analysis.sra_run:
        log.debug("transfering input file from sra")
        self.analysis.sra_run.download()
    else:
        log.debug("transfering input file from user")
        save_inputfile(analysis, inputfile)
    log.debug("done")
    return analysis

def save_inputfile(analysis, remote_file):
    local_file = open(analysis.inputfile_path, 'wb')
    from shutil import copyfileobj as copy
    copy(remote_file, local_file)
    local_file.close()

def start_analysis(analysis):
    url = "http://127.0.0.1:6543/analyses/%d"
    rq = urllib2.Request(url % analysis.id)
    rq.get_method = lambda: 'START'
    rsp = urllib2.urlopen(rq)
    body = rsp.read()
    rsp.close()
    if rsp.getcode() != 200:
        raise Exception("Worker failed to start analysis: " + body)

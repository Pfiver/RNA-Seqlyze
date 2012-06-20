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
    # create db object
    analysis = Analysis(**attributes)
    session.add(analysis)
    session.flush() # sets analysis.id
    analysis.create_directories()
    # transfer inputfile
    save_inputfile(analysis, inputfile)
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

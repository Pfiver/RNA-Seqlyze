import logging
log = logging.getLogger(__name__)

import os
import urllib2

import rnaseqlyze
from rnaseqlyze.core import security
from rnaseqlyze.core.orm import Analysis, User, RNASeqRun, UploadSession

def get_upload_session(db_session, id=None):
    if id:
        sess = db_session.query(UploadSession).get(id)
    else:
        sess = UploadSession()
        db_session.add(sess)
        db_session.flush()
        log.info("created upload session %d" % sess.id)
    return sess

def get_uploadfile(db_session, session, name, type):
    if not session.analysis:
        session.analysis = Analysis()
        db_session.add(session.analysis) # needed ?
        db_session.flush() # sets analysis.id
        session.analysis.create_data_dir()

    if type == 'inputfile':     # Short Reads in SRA or FASTQ format
        if session.analysis.inputfile_name:
            pass # TODO: remove
        session.analysis.inputfile_name = name
        log.debug("creating upload file '%s' for analysis #%d" % (
                                         name,             session.analysis.id))
        return open(session.analysis.inputfile_path, "w+b")
    elif type == 'genbankfile': # Organism genbank file
        if session.analysis.genbankfile_name:
            pass # TODO: remove
        session.analysis.genbankfile_name = name
        log.debug("creating upload file '%s' for analysis #%d" % (
                                         name,             session.analysis.id))
        return open(session.analysis.genbankfile_path, "w+b")
    else:
        raise Exception()

def create_analysis(db_session, attributes):
    # owner handling
    if 'owner' not in attributes:
        owner = db_session.query(User).get("anonymous")
        if not owner:
            owner = User("anonymous")
            db_session.add(owner)
        attributes['owner'] = owner
    # srr handling
    if 'rnaseq_run' in attributes:
        rnaseq_run = db_session.query(RNASeqRun).get(attributes['rnaseq_run'])
        if not rnaseq_run:
            try:
                rnaseq_run = RNASeqRun(attributes['rnaseq_run'])
                rnaseq_run.create_directories()
                db_session.add(rnaseq_run)
                attributes['rnaseq_run'] = rnaseq_run
            except:
                # The RNASeqRun constructor checks if the SRAnnnnnn argument
                # and raises an exception unless it passes the checks
                # e.g. if field was left blank/at default value
                del attributes['rnaseq_run']
                pass # TODO: decide/document what to do

    # get analysis
    # it already exists in case the user has uploaded something
    analysis = db_session.query(UploadSession) \
                    .get(attributes['upload_session']).analysis
    del attributes['upload_session']
    log.info(attributes)
    if not analysis:
        # create db object
        analysis = Analysis(**attributes)
        db_session.add(analysis)
        db_session.flush() # sets analysis.id
        analysis.create_data_dir()
    else:
        for attr, value in attributes.items():
            setattr(analysis, attr, value)

    # if no input file has been uploaded
    if not analysis.inputfile_name:
        # an SRR identifier is needed
        if not analysis.rnaseq_run:
            raise Exception("Please upload an input file or specify an SRR id")
        log.debug("transfering input file from sra")
        analysis.rnaseq_run.download()
        log.debug("done")

    return analysis

def start_analysis(analysis):
    url = "http://127.0.0.1:6543/analyses/%d"
    rq = urllib2.Request(url % analysis.id)
    rq.get_method = lambda: 'START'
    rsp = urllib2.urlopen(rq)
    body = rsp.read()
    rsp.close()
    if rsp.getcode() != 200:
        raise Exception("Worker failed to start analysis: " + body)

import logging
log = logging.getLogger(__name__)

import urllib2

import rnaseqlyze
from rnaseqlyze.core.entities import Analysis, User, RNASeqRun, UploadSession

def get_upload_session(db_session):
    sess = UploadSession()
    db_session.add(sess)
    db_session.flush()
    return sess

def get_uploadfile(db_session, session, name, type):
    # This doesn't look right, but it works. The database needs to
    # be locked here to make sure that the first upload request that
    # comes in creates the analysis and the second uses the same analysis.
    # We need to lock the whole database and this seemingly useless statement
    # does just that. With SQLite. I have been asking on irc #sqlalchemy about
    # how to do it the right way, but I didn't get any useful reply. I have
    # checked the SQLAlchemy as well as the SQLite docs and tried various things
    # like DBSession.execute("BEGIN") and such things - nothing seems to work
    # - this is the only solution I have found.
    #
    # EDT: it could be as simple as increasing the sqlalchemy debug level, check
    #      what sql statements are executed and then DBSession.execute() those
    #
    session.analysis = session.analysis
    #
    if not session.analysis:
        session.analysis = Analysis()
        db_session.add(session.analysis) # needed ?
        db_session.flush() # sets analysis.id
        session.analysis.create_data_dir()

    assert type in ('inputfile', 'genbankfile')

    typename = type + '_name'

    # inputfile_name -> Short Reads in SRA or FASTQ format
    # genabnkfile_name -> Organism genbank file
    #=if session.analysis.inputfile_name:
    if getattr(session.analysis, type + '_uploaded'):
        # you land here if a user uploads
        # more than one file per type
        # this is not intended, BUT
        # these are the interwebs!
        pass # FIXE: remove old

    #=session.analysis.inputfile_name = name
    setattr(session.analysis, typename, name)

    log.debug("creating upload file '%s' for analysis #%d" % (
                                     name,             session.analysis.id))

    # this would be the place to throw in a wrapper
    # to track upload progress the old way, i.e.
    # with server callbacks...

    #inputfile_path is a @property
    return open(getattr(session.analysis, type + '_path'), "w+b")

def get_analysis(db_session, attributes):

    # owner handling
    if 'owner' not in attributes:
        owner = db_session.query(User).get("anonymous")
        if not owner:
            owner = User("anonymous")
            db_session.add(owner)
        attributes['owner'] = owner

    # srr handling
    rnaseq_run = None
    if 'rnaseq_run' in attributes \
       and 'inputfile_name' not in attributes:
        log.debug("rnaseq_run: %s" % attributes['rnaseq_run'])
        rnaseq_run = db_session.query(RNASeqRun).get(attributes['rnaseq_run'])
        if rnaseq_run:
            attributes['rnaseq_run'] = rnaseq_run
        else:
            try:
                log.debug("creating new RNASeqRun")
                rnaseq_run = RNASeqRun(srr=attributes['rnaseq_run'])
                attributes['rnaseq_run'] = rnaseq_run
                rnaseq_run.create_directories()
                db_session.add(rnaseq_run)
            except Exception, e:
                # The RNASeqRun constructor checks the SRRnnnnnn argument
                # and raises an exception unless it passes the checks
                # e.g. if field was left blank/at default value
                # TODO: decide/document what to do
                log.debug("failed: %r" % e)
        del attributes['rnaseq_run']

    upload_session = db_session.query(UploadSession) \
                            .get(attributes['upload_session'])
    del attributes['upload_session']
    if not upload_session:
        raise Exception('this session has expired -'
                        ' reload the "New Analysis" page to start a new one')

    # the analysis exist already if the user uploaded something
    if upload_session.analysis:
        analysis = upload_session.analysis
        for attr, value in attributes.items():
            setattr(analysis, attr, value)

    else:
        # create db object
        log.debug("creating new analysis: %s" % attributes)
        analysis = Analysis(**attributes)
        db_session.add(analysis)
        db_session.flush() # sets analysis.id
        analysis.create_data_dir()

    # allow no more uploads to this analysis
    db_session.delete(upload_session)

    # if no input file has been uploaded
    if not analysis.inputfile_name:
        # an SRR identifier is needed
        if rnaseq_run:
            analysis.rnaseq_run = rnaseq_run
        else:
            raise Exception("Please upload an input file or specify an SRR id")

    return analysis

def start_analysis(analysis):
    rq = RNASWorkerSTARTRequest(
            "http://127.0.0.1:%s/analyses/%d" % (
                rnaseqlyze.worker_port, analysis.id))
    opener = urllib2.build_opener(HTTRNASWorkerHandler())
    rsp = opener.open(rq)
    body = rsp.read()
    rsp.close()

class RNASWorkerSTARTRequest(urllib2.Request):
    def get_method(self):
        return 'START'

class HTTRNASWorkerHandler(urllib2.HTTPHandler):
    def http_error(self, req, fp, code, msg, hdrs):
        raise WorkerException(fp.read())
    http_error_400 = http_error
    http_error_500 = http_error

class WorkerException(Exception):
    def __init__(self, exc_body):
        self.exc_body = exc_body
    def __repr__(self):
        return "WorkerException()"
    def __str__(self):
        return self.exc_body

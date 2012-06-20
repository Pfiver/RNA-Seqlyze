import logging
log = logging.getLogger(__name__)

import os
import urllib2

import rnaseqlyze
from rnaseqlyze.core import security
from rnaseqlyze.core.orm import Analysis, User

# TODO: make these functions methods of a mixin class

class AnalysisMixins(object):
    @property
	def data_dir(self):
        return os.path.join(rnaseqlyze.analyses_path, str(self.id))

    @property
	def shared_data_dir(self):
        acc = self.org_accession
        return os.path.join(rnaseqlyze.shared_data_path, acc)

    def create_directories(self)
        if not os.path.isdir(self.data_dir):
            os.makedirs(data_dir)
        if not os.path.isdir(self.shared_data_dir):
            os.makedirs(self.shared_data_dir)

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
        return "".join(fq_file.readline() for i in range(4))

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
    session.flush() # set analysis.id

    # transfer inputfile
    save_inputfile(analysis, inputfile)

    return analysis

def save_inputfile(analysis, remote_file):
    local_path = get_inputfile_path(analysis)
    local_file = open(local_path, 'wb')
    from shutil import copyfileobj as copy
    copy(remote_file, local_file)
    local_file.close()

def start_analysis(analysis):
    url = "http://127.0.0.1:6543/analyses/%d" % analysis.id
    STARTRequest = type('START', urllib2.Request, None)
    STARTRequest.get_method = lambda self: 'START'
    rsp = rllib2.urlopen(STARTRequest(url))
    body = rsp.read()
    rsp.close()
    if rsp.getcode() != 200:
        raise Exception("Worker failed to start analysis: " + body)


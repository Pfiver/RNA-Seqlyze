import os
import urllib2

import rnaseqlyze
from .orm import Analysis, User

datapath = rnaseqlyze.config.get("cache", "path")

class ServiceError(Exception):
    def __init__(self, msg, e):
        t = type(e)
        ServiceError.__init__(self,
            "%s (%s.%s(%s))" % (msg, t.__module__, t.__name__, e.args[0]))

def create_analysis(session, inputfile, attributes):

    # owner handling
    ################
    if 'owner' not in attributes:
        owner = session.query(User).get("anonymous")
        if not owner:
            owner = User("anonymous")
            session.add(owner)
        attributes['owner'] = owner

    # create db object
    ##################
    analysis = Analysis(**attributes)
    session.add(analysis)
    session.flush() # set analysis.id

    # transfer inputfile
    ####################
    save_inputfile(analysis, inputfile)

    return analysis


def save_inputfile(analysis, file):

    topdir = datapath + os.sep + 'analyses'
    if not os.path.isdir(topdir):
        os.mkdir(topdir)
    topdir += os.sep + str(analysis.id)
    if not os.path.isdir(topdir):
        os.mkdir(topdir)
    local_file = open(topdir + os.sep + "inputfile", 'wb')
    file.seek(0)
    while 1:
        data = file.read(4096)
        if not data: break
        local_file.write(data)
    local_file.close()


def start_analysis(analysis):
   try:
        url = "http://127.0.0.1:6543/analyses/%d" % analysis.id
        assert urllib2.urlopen(STARTRequest(url)).getcode() == 200
   except Exception, e:
       raise ServiceError("failed to notify worker", e)

class STARTRequest(urllib2.Request):
    def get_method(self):
        return 'START'

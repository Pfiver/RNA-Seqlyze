import os

import rnaseqlyze
from .orm import Analysis, User

datapath = rnaseqlyze.config.get("cache", "path")

class ServiceError(Exception):
    def __init__(self, msg, e):
        t = type(e)
        msg = "%s (%s.%s)" % (
            msg, t.__module__, t.__name__)
        super(ServiceError, self).__init__(msg)

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

        con = urllib2.urlopen(STARTRequest(url))
        con.read()
        con.close()

        assert response['status'] == 200

    except Exception, e:
        raise ServiceError("failed to notify worker", e)

import urllib2
class STARTRequest(urllib2.Request):
    def get_method(self):
        return 'START'

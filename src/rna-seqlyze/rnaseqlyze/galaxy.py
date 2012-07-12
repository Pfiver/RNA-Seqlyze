"""
RNA-Seqlyze Galaxy Module

Shamelessly piggy-back onto Penn-State University's "Galaxy" Project.

RNA-Seqlyze needs a some publicly available Web-Space, which PSU provides
plenty of for bioinformatics reseach data (250.0 Gb per user as of 4 Jul 2012).

Thanks go to Penn-State University!

http://www.psu.edu/

"""

# FIXME: the whole code here needs heavy refactoring

import logging
log = logging.getLogger(__name__)

import os, json, time, ftplib, \
       urllib, urllib2, cookielib
from time import time
from threading import local
from datetime import datetime, timedelta

import lxml.html

import rnaseqlyze
from rnaseqlyze import multipart

email = 'ucgxccgr@mailinator.com'
password = 'brtbhcdg'
api_key = 'dddb2c53c96c0c4d263e6c74b507d203'
hostname = 'main.g2.bx.psu.edu'

default_history = '16f9a8e916e0e908'

default_history_url = 'https://main.g2.bx.psu.edu/u/dcgdftvcdv/h/rna-seqlyze'

history_path_template = '/api/histories/{history}/contents'
ucsc_bam_track_template = \
        '/display_application/{dataset}/ucsc_bam/archaea/None/param/track'

ucsc_bam_path_template = \
        '/display_application/{dataset}/' \
        'ucsc_bam/archaea/None/data/galaxy_{dataset}.bam'

dataset_info_url_template = "/api/histories/{history}/contents/{dataset}"

dataset_display_url_template = "/datasets/{dataset}/display"

rq_headers = {}

class Session(local):
    cookies = None
    created = None

session = Session()

def api_call(path):
    url = "https://" + hostname + path
    return urllib2.urlopen(url + "?key=" + api_key).read()

def login():
    cookie_jar = cookielib.CookieJar()
    urllib2.install_opener(urllib2.build_opener(
                    urllib2.HTTPCookieProcessor(cookie_jar)))
    log.info("Loggin in to galaxy server %s ..." % hostname)
    login = "https://" + hostname + "/user/login"
    rq = urllib2.Request(login, headers=rq_headers)
    request = urllib2.urlopen(rq)
    doc = lxml.html.parse(request).getroot()
    form = doc.forms[0]
    form.fields["email"] = email
    form.fields["password"] = password
    submit = "login_button", form.fields["login_button"]
    data = urllib.urlencode(form.form_values() + [submit])
    log.debug("posting login form: %s" % form.action)
    rq = urllib2.Request(form.action, headers=rq_headers)
    request = urllib2.urlopen(rq, data)
    doc = lxml.html.parse(request).getroot()
    log.info("Success!")
    return cookie_jar

def import_upload(filename):
    if not (session.created
            and session.created > (datetime.now() - timedelta(minutes=30))):
        session.cookies = login()
        session.created = datetime.now()

    urllib2.install_opener(urllib2.build_opener(
                    urllib2.HTTPCookieProcessor(session.cookies)))
    log.info("Importing ftp file")
    tool = "https://" + hostname + "/tool_runner?tool_id=upload1"
    rq = urllib2.Request(tool, headers=rq_headers)
    request = urllib2.urlopen(rq)
    doc = lxml.html.parse(request).getroot()
    found = False
    form = doc.forms[0]
    inp = form.inputs["files_0|ftp_files"]
    if isinstance(inp, lxml.html.InputElement):
        if inp.attrib['value'] == filename:
            found = inp.checked = True
    elif isinstance(inp, lxml.html.CheckboxGroup):
        for box in inp:
            if box.attrib['value'] == filename:
                found = box.checked = True
    else:
        raise Exception("unexpected html element: %s" % inp)
    if not found:
        raise Exception("file not available for import: %s" % filename)
    submit = "runtool_btn", form.fields["runtool_btn"]
    data = multipart.urlencode(form.form_values() + [submit])
    log.debug("posting upload form: %s" % form.action)
    rq = urllib2.Request(form.action, data, headers=rq_headers)
    request = multipart.urlopen(rq)
    doc = lxml.html.parse(request).getroot()
    log.info("Success!")

def ftpupload(fileobj, filename):
    """
    upload a file object to galaxy
    based on http://love-python.blogspot.com/2008/02/ftp-file-upload.html
    """
    log.info("uploading file to ftp server")
    ftp = ftplib.FTP(hostname, email, password)

    try:
        total = os.stat(fileobj.name).st_size
    except:
        total = 0

    class cbc(object):
        def __init__(self, total):
            self.total = total
            self.sent = 0
            self.then = time()
        def __call__(self, buf):
            self.sent += len(buf)
            now = time()
            if self.then < now - 15:
                self.then = now
                log.info("%d kb left" % ((self.total - self.sent) / 1024))

    log.info("sending %d kb of data..." % ((total / 1024 or -1)))
    ftp.storbinary('STOR ' + filename, fileobj, callback=cbc(total))
    log.info("Success!")
    ftp.quit()

def upload(fileobj, filename):

    # can't initialize this at module import time
    # because rnaseqlyze.xxx properties not initialized
    global rq_headers
    try:
        mail = rnaseqlyze.admin_email
    except:
        # rnaseqlyze not .configure()d
        mail = os.getenv("USER") + "@" + os.uname()[1]
    rq_headers = {
        'User-Agent': "%s (version:%s / admin:%s)" % (
            rnaseqlyze.project_name, rnaseqlyze.__version__, mail),
    #    'User-Agent': "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:13.0)" \
    #                  " Gecko/20100101 Firefox/13.0.1",
    }

    ftpupload(fileobj, filename)
    import_upload(filename)
    datasets = json.loads(api_call(
        history_path_template.format(history=default_history)))
    # assume objects are ordered chronologically...
    for dataset in reversed(datasets):
        if dataset['name'] == filename:
            return dataset['id']
    raise Exception("Couldn't find id of uploaded file in dataset")

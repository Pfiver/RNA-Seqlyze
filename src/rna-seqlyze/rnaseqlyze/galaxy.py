import logging
log = logging.getLogger(__name__)

import os
import urllib
import urllib2
import lxml.html
import cookielib
import multipart

email = 'ucgxccgr@mailinator.com'
password = 'brtbhcdg'
api_key = 'dddb2c53c96c0c4d263e6c74b507d203'
hostname = 'main.g2.bx.psu.edu'

default_history = 'b8468b3367a258a6'
default_history_url = 'https://main.g2.bx.psu.edu/u/dcgdftvcdv/h/biocalc'

history_path_template = '/api/histories/%(history)s/contents'
ucsc_bam_track_template = \
        '/display_application/%(dataset)s/ucsc_bam/archaea/None/param/track'

def api_call(path):
    url = "https://" + hostname + path
    return urllib2.urlopen(url + "?key=" + api_key).read()

def login(cookie_file=None):
    if not cookie_file:
        cookie_jar = cookielib.CookieJar()
    else:
        cookie_jar = cookielib.MozillaCookieJar()
    urllib2.install_opener(urllib2.build_opener(
                                    urllib2.HTTPCookieProcessor(cookie_jar)))
    log.info("Loggin in to galaxy server %s ..." % hostname)
    login = "https://" + hostname + "/user/login"
    request = urllib2.urlopen(login)
    doc = lxml.html.parse(request).getroot()
    form = doc.forms[0]
    form.fields["email"] = email
    form.fields["password"] = password
    submit = "login_button", form.fields["login_button"]
    data = urllib.urlencode(form.form_values() + [submit])
    log.debug("posting login form: %s" % form.action)
    request = urllib2.urlopen(form.action, data)
    doc = lxml.html.parse(request).getroot()
    log.info("Success!")
    if cookie_file:
        cookie_jar.save(cookie_file, ignore_discard=True, ignore_expires=True)
    return cookie_jar

def import_uploads(cookie_jar=None, cookie_file="cookies.txt"):
    if cookie_jar:
        cookie_file = None
    else:
        cookie_jar = cookielib.MozillaCookieJar()
        cookie_jar.load(cookie_file, ignore_discard=True, ignore_expires=True)
    urllib2.install_opener(urllib2.build_opener(
                                    urllib2.HTTPCookieProcessor(cookie_jar)))
    log.info("Importing ftp files...")
    tool = "https://" + hostname + "/tool_runner?tool_id=upload1"
    request = urllib2.urlopen(tool)
    doc = lxml.html.parse(request).getroot()
    form = doc.forms[0]
    inp = form.inputs["files_0|ftp_files"]
    if isinstance(inp, lxml.html.InputElement):
        inp.checked = True
    elif isinstance(inp, lxml.html.CheckboxGroup):
        for group in inp:
            group.checked = True
    else:
        raise Exception("unexpected html element: %s" % inp)
    submit = "runtool_btn", form.fields["runtool_btn"]
    data = multipart.urlencode(form.form_values() + [submit])
    log.debug("posting upload form: %s" % form.action)
    request = multipart.urlopen(form.action, data)
    doc = lxml.html.parse(request).getroot()
    log.info("Success!")
    if cookie_file:
        cookie_jar.save(cookie_file, ignore_discard=True, ignore_expires=True)

def ftpupload(fileobj, filename):
    """
    upload a file object to galaxy
    based on http://love-python.blogspot.com/2008/02/ftp-file-upload.html
    """
    import os.path, ftplib
    log.debug("uploading file to ftp server")
    ftp = ftplib.FTP(hostname, email, password)
    ftp.storbinary('STOR ' + filename, fileobj)
    log.debug("Success!")
    ftp.quit()

def upload(fileobj, filename):
    import json
    ftpupload(fileobj, filename)
    import_uploads(login())
    histories = json.loads(api_call(
        history_path_template % dict(history=default_history)))
    # assume objects are ordered chronologically...
    for history in reversed(histories):
        if history['name'] == filename:
            return history['id']

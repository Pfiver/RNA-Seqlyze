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

default_history = '6556fe2755b424cd'
history_path_template = '/api/histories/$history/contents'
ucsc_bam_track_template = '/display_application/$dataset/ucsc_bam/archaea/$history/param/track'

def api_call(path):
    url = "https://" + hostname + path
    return urllib2.urlopen(url + "?key=" + api_key).read()

def login():
    cookie_file = "cookies.txt"
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
    cookie_jar.save(cookie_file, ignore_discard=True, ignore_expires=True)

def import_uploads():
    cookie_file = "cookies.txt"
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
    cookie_jar.save(cookie_file, ignore_discard=True, ignore_expires=True)

def upload(fileobj):
    """
    upload a file object to galaxy
    based on http://love-python.blogspot.com/2008/02/ftp-file-upload.html
    """
    import os.path, ftplib
    sftp = ftplib.FTP(hostname, email, password)
    sftp.storbinary('STOR ' + os.path.basename(file), fileobj)
    sftp.quit()

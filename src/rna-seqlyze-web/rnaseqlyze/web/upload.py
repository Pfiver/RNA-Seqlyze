"""
Pyramid Application Upload View

This module handles the upload of analysis files.

The upload interface consists of `plupload`
(from http://www.plupload.com/) on the client
and this hack on the server side. Combining the
two and creating a working solution was not trivial.

Documentation and inspiration to create this
was, amongst others, taken from the following documents:

 - http://www.plupload.com/documentation.php
 - http://hg.python.org/cpython/file/2.7/Lib/cgi.py#l353
 - https://raw.github.com/moxiecode/plupload/master/examples/upload.php
 - https://github.com/hcwebdev/plupload/blob/master/examples/server.py
 - https://github.com/Pylons/webob/blob/master/webob/request.py#L102
 - https://hg.gawel.org/gp.fileupload/file/default/gp/fileupload/storage.py#l97
"""

import logging
log = logging.getLogger(__name__)

import transaction
from pyramid.view import view_config

from rnaseqlyze.web import DBSession
from rnaseqlyze.core import service
from rnaseqlyze.core.orm import UploadSession

@view_config(route_name='upload', request_method='POST', renderer="json")
def upload(request):
    log.debug("upload(): content-type '%s'" % request.content_type)
    fs = FieldStoragx(fp=request.environ['wsgi.input'], environ=request.environ)
    return dict(jsonrpc="2.0", result=None, id=None)

import cgi
class FieldStoragx(cgi.FieldStorage):
    def __init__(self, fp=None, headers=None, outerboundary="",
                 environ=None, keep_blank_values=0, strict_parsing=0):
        self.environ = environ
        cgi.FieldStorage.__init__(self, fp, headers, outerboundary,
                                  environ, keep_blank_values, strict_parsing)
        if self.filename:
            return
        assert len(self.value) < 1000
        if self.name == 'session':
            environ['rnaseqlyse.upload_session'] = \
                   DBSession.query(UploadSession).get(int(self.value))
        elif self.name in ('name', 'type'):
            environ['rnaseqlyse.upload_' + self.name] = self.value
        else:
            return
        log.debug("FieldStoragx(%s -> %s)" % (self.name, self.value))

    def make_file(self, binary=None):

        assert self.filename
        log.debug("FieldStoragx.make_file(%s)" % self.filename)

        args = {}
        for kw in 'session', 'name', 'type':
            args[kw] = self.environ['rnaseqlyse.upload_' + kw]

        fd = service.get_uploadfile(DBSession, **args)
        # commit the (managed) session early here, so later
        # requests can re-use the Analysis object that the
        # first one has implicitly created by calling
        # service.get_uploadfile
        import transaction
        transaction.commit()

        return fd


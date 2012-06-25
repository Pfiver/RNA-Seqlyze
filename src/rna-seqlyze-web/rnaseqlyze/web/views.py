import logging
log = logging.getLogger(__name__)

from pyramid.view import view_config
from pyramid.response import Response
from pyramid.renderers import get_renderer
from pyramid.httpexceptions import (
        HTTPFound, HTTPError, HTTPServiceUnavailable, HTTPInternalServerError
)

import transaction
from sqlalchemy.exc import DBAPIError

from rnaseqlyze.web import DBSession
from rnaseqlyze.core import service
from rnaseqlyze.core.orm import Analysis

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
                   service.get_upload_session(DBSession, self.value)
        elif self.name in ('name', 'type'):
            environ['rnaseqlyse.upload_' + self.name] = self.value
        else:
            return
        log.debug("FieldStoragx(%s -> %s)" % (self.name, self.value))

    def make_file(self, binary=None):
        args = {}
        assert self.filename
        for kw in 'session', 'name', 'type':
            args[kw] = self.environ['rnaseqlyse.upload_' + kw]
        log.debug("FieldStoragx.make_file(%s)" % self.filename)
        fd = service.get_uploadfile(DBSession, **args)
        import transaction
        transaction.commit()
        return fd

@view_config(route_name='analyses', renderer='templates/create.pt')
def create(request):
    sess = service.get_upload_session(DBSession)
    return { 'upload_session': sess.id }

@view_config(route_name='analyses', request_method='POST')
def post(request):

    # TODO: csrf security checks
    #       see "shootout" pyramid demo app
    try:
        analysis = service.get_analysis(
                    DBSession.unmanaged, attributes=request.POST)

        # the analysis must exist in the database
        # so the worker can find it and start working
        DBSession.unmanaged.commit()

        service.start_analysis(analysis)
        log.debug("started analysis #%d by '%s'" % (
                            analysis.id, analysis.owner.name))

        return HTTPFound(request.route_path('analysis', id=analysis.id))
    except:
        log.info("abort")
        DBSession.unmanaged.rollback()
        log.debug("rollback complete")
        raise

@view_config(route_name='analysis', renderer='templates/analysis.pt')
def display(request):
    id = int(request.matchdict["id"])
    return { 'analysis': DBSession.query(Analysis).get(id) }

@view_config(context=DBAPIError)
def dbapi_error(request):
    import traceback
    detail = conn_err_msg
    detail += traceback.format_exc(999)
    log.debug(detail)
    return HTTPInternalServerError(detail=detail)

@view_config(context=Exception)
def error(request):
    import traceback
    detail = traceback.format_exc(999)
    log.debug(detail)
    return HTTPInternalServerError(detail=detail)

import string
HTTPError.body_template_obj = string.Template(
    "${explanation}${br}${br}\n<pre>${detail}</pre>\n${html_comment}\n\n")

conn_err_msg = """\
Pyramid is having
a problem connecting to to database.

You may need to run the
"initialize_rnaseqlyze_db" script to create it.

Afterwards, restart the Pyramid application, i.e. send a
SIG_INT to the apache mod_wsgi daemon processes, and try again.

"""

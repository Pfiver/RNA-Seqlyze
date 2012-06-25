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
    log.debug(request.content_type)
    fs = FieldStoragx.from_request(request)
    return dict(jsonrpc="2.0", result=None, id=None)

import cgi
class FieldStoragx(cgi.FieldStorage):
    @classmethod
    def from_request(request):
        wsgi_input = request.environ['wsgi.input']
        return FieldStoragx(fp=wsgi_input, environ=request.environ)

    def __init__(self, *args, **kwargs):
        cgi.FieldStorage.__init_(self, *args, **kwargs)
        log.debug("FieldStoragx: %s -> %s" % (self.name, self.value))

        if self.filename:
            args = {}
            for kw in 'session', 'name', 'type':
                args[kw] = self.environ['rnaseqlyse.upload_' + info]
            self.uploadfile = service.get_uploadfile(DBSession, **args)

        assert len(self.value) < 1000
        if self.name = 'upload_session':
            self.environ['rnaseqlyse.upload_session'] = \
                          service.get_upload_session(DBSession, self.value)
        elif self.name in ('name', 'type'):
            self.environ['rnaseqlyse.upload_' + self.name] = self.value
        else
            raise

    def make_file(self, binary=None):
        if not hasattr(self, 'uploadfile'):
            raise "Unexpected input"
        return self.uploadfile

@view_config(route_name='analyses', renderer='templates/create.pt')
def create(request):
    import sha, datetime
    return { 'session': sha.new(str(datetime.datetime.now())).hexdigest() }

@view_config(route_name='analyses', request_method='POST')
def post(request):

    # TODO: csrf security checks
    #       see "shootout" pyramid demo app
    post = request.POST

    # inputfile handling
    try:
        # .filename attribute access will fail if no file is specified
        post['inputfile_name'] = post['inputfile'].filename
        inputfile = post['inputfile'].file
    except:
        inputfile = None
        # the user should have provided a
        # meaningful SRR identifier in sra_run in this case
        pass
    del post['inputfile']

    # organism handling
    if 'org_accession' not in post:
        if 'genebankfile' in post:
            raise Exception('Genebank file upload not yet implemented')
        else:
            raise Exception('Organism not specified')
    else:
        if 'genebankfile' in post:
            del post['genebankfile']

    analysis = service.create_analysis(DBSession.unmanaged,
                                       inputfile, attributes=post)
    DBSession.unmanaged.commit()
    # the analysis must exist in the database
    # so the worker can find it and start working

    service.start_analysis(analysis)
    log.debug("started analysis #%d by '%s'" % (
                        analysis.id, analysis.owner.name))

    return HTTPFound(request.route_path('analysis', id=analysis.id))

@view_config(route_name='analysis', renderer='templates/analysis.pt')
def display(request):
    id = int(request.matchdict["id"])
    return { 'analysis': DBSession.query(Analysis).get(id) }

@view_config(context=DBAPIError)
def dbapi_error(request):
    import traceback
    detail = conn_err_msg
    detail += traceback.format_exc(999)
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

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

@view_config(route_name='analyses', request_method='POST')
def post(request):

    # TODO: csrf security checks
    #       see "shootout" pyramid demo app
    post = request.POST

    # inputfile handling
    post['inputfile_name'] = post['inputfile'].filename
    inputfile = post['inputfile'].file
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
    # for the worker to be able to start working on it

    service.start_analysis(analysis)
    log.debug("started analysis #%d by '%s'" % (
                        analysis.id, analysis.owner.name))

    return HTTPFound(request.route_path('analysis', id=analysis.id))


@view_config(route_name='analysis', renderer='templates/analysis.pt')
def display(request):
    id = int(request.matchdict["id"])
    return {
        'service': service,
        'analysis': DBSession.query(Analysis).get(id),
    }

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

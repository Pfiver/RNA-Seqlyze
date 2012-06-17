from pyramid.view import view_config, view_defaults
from pyramid.response import Response
from pyramid.renderers import get_renderer
from pyramid.httpexceptions import HTTPFound

from sqlalchemy.exc import DBAPIError

from . import DBSession
from ..core import service
from ..core.orm import Analysis

import logging
log = logging.getLogger(__name__)


@view_config(route_name='home', renderer='templates/home.pt')
def home(request):
    return {}


@view_config(route_name='analyses', renderer='templates/create.pt')
def create(request):
    return {}


@view_config(route_name='analyses', request_method='POST')
def post(request):

    # TODO: csrf security checks
    #       see "shootout" pyramid demo app
    try:
        post = request.POST

        # inputfile handling
        ####################
        post['inputfilename'] = post['inputfile'].filename
        inputfile = post['inputfile'].file
        del post['inputfile']

        # organism handling
        ###################
        if 'org_accession' not in post:
            if 'genebankfile' in post:
                raise Exception('Genebank file upload not yet implemented')
            else:
                raise Exception('Organism not specified')
        else:
            if 'genebankfile' in post:
                del post['genebankfile']

        analysis = service.create_analysis(DBSession, inputfile, attributes=post)

        log.debug("starting analysis #%d by '%s'..." % (analysis.id, analysis.owner.name))

        service.start_analysis(analysis)

    except DBAPIError, e:
        return Response(conn_err_msg % e, content_type='text/plain', status_int=500)

    return HTTPFound(request.route_path('analysis', id=str(analysis.id)))


@view_config(route_name='analysis', renderer='templates/analysis.pt')
def display(request):

    try:
        id = int(request.matchdict["id"])
        analysis = DBSession.query(Analysis).get(id)

    except DBAPIError, e:
        return Response(conn_err_msg % e, content_type='text/plain', status_int=500)

    return {
        'analysis': analysis,
    }


conn_err_msg = """\
Pyramid is having
a problem connecting to to database.

You may need to run the
"initialize_rnaseqlyze_db" script to create it.

Afterwards, restart the Pyramid application, i.e. send a
SIG_INT to the apache mod_wsgi daemon processes, and try again.

The error was: %s
"""

from pyramid.view import view_config
from pyramid.response import Response

from sqlalchemy.exc import DBAPIError

from rnaseqlyze.core import service
from rnaseqlyze.core.orm import Analysis
from rnaseqlyze.web import DBSession


@view_config(route_name='home', renderer='templates/home.pt')
def home(request):
    return {}


from pyramid.httpexceptions import HTTPFound
@view_config(route_name='new')
def new(request):
    try:
        analysis = service.create_anonymous_analysis(DBSession)
        DBSession.add(analysis)
        DBSession.flush() # set analysis.id
        loc = request.route_path('analysis', id=str(analysis.id))
        return HTTPFound(location=loc)
    except DBAPIError, e:
        return Response(conn_err_msg % e, content_type='text/plain', status_int=500)


@view_config(route_name='analysis', renderer='templates/analysis.pt')
def analysis(request):
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

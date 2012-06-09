from pyramid.view import view_config
from pyramid.response import Response

from sqlalchemy.exc import DBAPIError

from rnaseqlyze.core.orm import DBSession, Analysis


@view_config(route_name='home', renderer='templates/home.pt')
def home(request):
    return {}


@view_config(route_name='new', renderer='templates/analysis.pt')
@view_config(route_name='analysis', renderer='templates/analysis.pt')
def analysis(request):
    try:
        analyses_query = DBSession.query(Analysis)
        analyses = analyses_query.all()
    except DBAPIError:
        return Response(conn_err_msg, content_type='text/plain', status_int=500)

    return {
        'analyses': analyses,
    }


conn_err_msg = """\
Pyramid is having
a problem connecting to to database.

You may need to run the
"initialize_rnaseqlyze_db" script to create it.

Afterwards, restart the Pyramid application, i.e. send a
SIG_INT to the apache mod_wsgi daemon processes, and try again.
"""

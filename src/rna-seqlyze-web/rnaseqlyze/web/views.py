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


class Analysis(object):

    def __init__(self, request):
        self.request = request
        renderer = get_renderer("templates/analysis-base.pt")
        self.analysis_base = renderer.implementation()


    @view_config(route_name='analysis', request_method='GET', renderer='templates/analysis-input.pt')
    def input(request):
        try:
            id = int(request.matchdict["id"])
            analysis = DBSession.query(Analysis).get(id)
        except DBAPIError, e:
            return Response(conn_err_msg % e, content_type='text/plain', status_int=500)

        return {
            'analysis': analysis,
            'pagescript': 'analysis.js',
        }


    @view_config(route_name='analysis', request_method='POST', renderer='templates/analysis-process.pt')
    def process(request):
        try:
            id = int(request.matchdict["id"])
            analysis = DBSession.query(Analysis).get(id)
            input_file = request.POST["inputfile"].file
            analysis.inputfilename = request.POST["inputfile"].filename
            file_path = os.path.join('/home/pfeifer/data/cache/analyses/%d/%s' % (id, "inputfile"))
            output_file = open(file_path, 'wb')
            input_file.seek(0)
            while 1:
                data = input_file.read(2<<12)
                if not data:
                    break
                output_file.write(data)
            output_file.close()

        except DBAPIError, e:
            return Response(conn_err_msg % e, content_type='text/plain', status_int=500)

        return {
            'analysis': analysis,
            'pagescript': 'analysis.js',
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

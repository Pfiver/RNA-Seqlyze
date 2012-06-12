from pyramid.view import view_config, view_defaults
from pyramid.response import Response
from pyramid.renderers import get_renderer

from sqlalchemy.exc import DBAPIError

from . import DBSession, datapath
from ..core import service
from ..core.orm import Analysis


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


@view_defaults(route_name='analysis')
class Analyses(object): # call the view here Analyses to avoid
                        # name collision with domain object Analysis

    def __init__(self, request):
        self.request = request
        self.base = get_renderer("templates/analysis-base.pt").implementation()


    def finished(context, request):
        id = int(request.matchdict["id"])
        return DBSession.query(Analysis).get(id).started


    @view_config(request_method='GET',
                 renderer='templates/analysis-input.pt')
    def input(self):
        try:
            id = int(self.request.matchdict["id"])
            analysis = DBSession.query(Analysis).get(id)
        except DBAPIError, e:
            return Response(conn_err_msg % e, content_type='text/plain', status_int=500)

        return {
            'analysis': analysis,
        }


    @view_config(request_method='GET',
                 custom_predicates=[finished],
                 renderer='templates/analysis-results.pt')
    def results(self):
        try:
            id = int(self.request.matchdict["id"])
            analysis = DBSession.query(Analysis).get(id)
        except DBAPIError, e:
            return Response(conn_err_msg % e, content_type='text/plain', status_int=500)

        return {
            'analysis': analysis,
        }


    @view_config(request_method='POST', renderer='templates/analysis-process.pt')
    def process(self):
        try:
            import os
            id = int(self.request.matchdict["id"])
            topdir = datapath + os.sep + 'analyses'
            analysis = DBSession.query(Analysis).get(id)
            if analysis.started:
                raise Exception("analysis already started")
            post = self.request.POST
            input_file = post["inputfile"].file
            analysis.inputfilename = post["inputfile"].filename
            analysis.organism = post["organism"]
            analysis.strandspecific = "strandspecific" in post
            analysis.pairended = "pairended" in post
            analysis.pairendlen = post["pairendlen"]
            if not os.path.isdir(topdir):
                os.mkdir(topdir)
            topdir += os.sep + str(id)
            if not os.path.isdir(topdir):
                os.mkdir(topdir)
            output_file = open(topdir + os.sep + "inputfile", 'wb')
            input_file.seek(0)
            while 1:
                data = input_file.read(2<<12)
                if not data:
                    break
                output_file.write(data)
            output_file.close()

            analysis.started = True

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

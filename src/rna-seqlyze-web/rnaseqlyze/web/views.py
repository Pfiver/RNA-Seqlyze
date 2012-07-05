"""
Pyramid Application User Views
"""

import logging
log = logging.getLogger(__name__)

from os.path import join

from pyramid.view import view_config
from pyramid.response import FileResponse
from pyramid.httpexceptions import HTTPFound

import transaction
from sqlalchemy.exc import DBAPIError

import rnaseqlyze
from rnaseqlyze.web import DBSession, DBSession_unmanaged
from rnaseqlyze.core import service
from rnaseqlyze.core.entities import Analysis, UCSCOrganism

@view_config(route_name='home', renderer='templates/home.pt')
def home(request):
    """
    **Home Page**

    This is the main entry point to the application. I.e. the landing page,
    the page that users see first.
    """
    return {}

@view_config(route_name='analyses', renderer='templates/create.pt')
def create(request):
    """
    **Create Page**

    This page is shown when the "New Analysis" button is clicked.
    """
    sess = service.get_upload_session(DBSession)
    return { 'upload_session': sess.id }

@view_config(route_name='analysis', renderer='templates/analysis.pt')
def display(request):
    """
    **Analysis Page**

    This page is displayed after the the anaysis has been created.
    When the user clicks "Submit" on the 'create' page, after
    the files are uploaded and the form information is submitted
    to the :func:`~post` view, the browser is redirected here.

    The page can also be viewed any time later on, no matter
    weather the analysis has already been completed or not.

    In case it is not yet completed, the page is constantly
    updated via XMLHTTPRequests to reflect the current status.
    """
    
    id = int(request.matchdict["id"])
    return { 'analysis': DBSession.query(Analysis).get(id) }

@view_config(route_name='analysis_files')
def analysis_files(request):
    """
    **Files View**

    This view serves up the files associated with
    an analysis on 'http://<rnaseqlyze>/analysis/{id}/files'.
    """
    return FileResponse(join(rnaseqlyze.analyses_path,
                request.matchdict['id'], *request.subpath))

import mimetypes
#mimetypes.add_type("text/plain", ".")
#mimetypes.add_type("text/plain", ".")
#mimetypes.add_type("text/plain", ".")
#mimetypes.add_type("text/plain", ".")
mimetypes.add_type("text/plain", ".gb")
mimetypes.add_type("text/plain", ".log")
mimetypes.add_type("text/plain", ".log0")
mimetypes.add_type("text/plain", ".info")
# FileResponse automatically sets the Content-Type header based on this

@view_config(route_name='analyses', request_method='POST')
def post(request):
    """
    **Create-Form Action**

    This view just redirects the client to the created analysis page.
    Before it is actually called, the files to be analyzed, are uploaded
    using the :func:`~.upload.upload` view callable.
    """
    # for documentation on the documentation reference syntax, see
    # http://sphinx.pocoo.org/domains.html#cross-referencing-python-objects

    # TODO: csrf security checks
    #       see "shootout" pyramid demo app

    # note:
    #  when using the "DBSession" (managed), the
    #  try:/except: rollback construct is not needed
    #  because the session is automatically rolled back
    #  otoh, if the _unmanaged session is used, it _has_ to
    #  be manually commited or rolled back if objects are modified

    try:
        analysis = service.get_analysis(
                    DBSession_unmanaged, attributes=request.POST)

        # the analysis must exist in the database
        # so the worker can find it and start working
        DBSession_unmanaged.commit()

        service.start_analysis(analysis)
        log.debug("started analysis #%d by '%s'" % (
                            analysis.id, analysis.owner.name))

        return HTTPFound(request.route_path('analysis', id=analysis.id))
    except:
        log.info("abort")
        DBSession_unmanaged.rollback()
        log.debug("rollback complete")
        raise

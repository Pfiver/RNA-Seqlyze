"""
Pyramid Application Views

To understand the applicatios architecture, head over to the
wonderful world of the pyramid web framework at http://www.pylonsproject.org/.

In case you have trouble with anything pyramid-related,
use the `source code on github <https://github.com/Pylons/pyramid>`_
or ask 'mcdonc' on freenode irc
`#pyramid <http://webchat.freenode.net/?channels=#pyramid>`_.
"""

import logging
log = logging.getLogger(__name__)

from pyramid.view import view_config
from pyramid.response import Response, FileResponse
from pyramid.httpexceptions import (
        HTTPFound, HTTPError, HTTPServiceUnavailable, HTTPInternalServerError
)

import transaction
from sqlalchemy.exc import DBAPIError

from rnaseqlyze.web import DBSession
from rnaseqlyze.core import service
from rnaseqlyze.core.orm import Analysis

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
    import os
    return FileResponse(
            os.path.join(rnaseqlyze.analyses_path,
                request.matchdict['id'], *request.subpath))

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

    # when using the "DBSession" (managed), the
    # try:/except: rollback construct is not needed
    # because the session is automatically rolled back

    # otoh, if the .unmanaged session is used, it _has_ to
    # be manually commited or rolled back if objects are modified
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

@view_config(context=Exception)
def error(request):
    """
    **Exception view**

    This is a catch-all view that serves up any errors
    that have occured while processing the a request.
    """
    detail = "An Exception was raised in rnaseqlyze.web: %s" % request.context

    if log.getEffectiveLevel() <= log.DEBUG:
        if isinstance(request.context, DBAPIError):
            detail += """
                This is a database related eror.

                If it is not yet initialized or the schema has changed,
                just run the "rnas-dbinit" script to (re-)initialize it.

                Afterwards, restart the Pyramid application, i.e. send a
                SIG_INT to the apache mod_wsgi daemon processes, and try again.

                """
        import traceback
        detail += '\n' + traceback.format_exc(999)

    log.error(detail)
    return HTTPInternalServerError(detail=detail)

# `Patch` HTTPError to get custom error messages (adding the <pre> tags)
from string import Template
HTTPError.body_template_obj = Template(
    "${explanation}${br}${br}\n<pre>${detail}</pre>\n${html_comment}\n\n")

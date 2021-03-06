"""
Pyramid Application Custom Error Views
"""

import logging
log = logging.getLogger(__name__)

from string import Template

from pyramid.view import view_config
from pyramid.response import Response, FileResponse
from pyramid.httpexceptions import (
        HTTPFound, HTTPError, HTTPServiceUnavailable, HTTPInternalServerError
)

import transaction
from sqlalchemy.exc import DBAPIError

import rnaseqlyze
from rnaseqlyze.web import DBSession, DBSession_unmanaged
from rnaseqlyze.core import service
from rnaseqlyze.core.entities import Analysis

@view_config(context=Exception)
def error(request):
    """
    **Exception view**

    This is a catch-all view that serves up any errors
    that have occured while processing the a request.

    The view just creates and returns a custom error response object.
    """
    return HTTPRNASeqError(request.exc_info)

class HTTPRNASeqError(HTTPError):
    """
    Custom HTTP Error class.

    This is a custom HTTP error class that extends
    :class:`pyramid.httpexceptions.HTTPError`, which extends
    :class:`pyramid.httpexceptions.WSGIHTTPException`. Have a look at the
    `source code <http://git.io/CqrfOg#L157>`_ to see how it works.

    It's ``code`` is 500, which generaly means "Internal Server Error".  If the
    application is in debugging mode -- i.e. the log level is DEBUG or less, a
    stack trace is added to the generated page as well as the log file.
    Otherwise, an informational message is displayed and only one line,
    containing the type of the error is logged.
    """
    code = 500
    title = "RNA-Seqlyze Web Application Error"
    explanation = "An Exception was raised in rnaseqlyze.web"
    html_template_obj = Template(Template('\n'.join(map(lambda s: s[8:], """\
        <html>
        <head>
        <title>${title}</title>
        </head>
        <body style="margin: 20px;">
        <h1>${title}</h1>
        ${body}
        </body>
        </html>
        """.split('\n')))).safe_substitute(title=title))
    def __init__(self, exc_info):
        e = exc_info[1]
        log.error(repr(e))
        body_template = "<b>${explanation}</b>\n<hr/>\n"
        cls = e.__class__.__name__
        if not e.args:
            self.explanation = "%s" % cls
        else:
            self.explanation = "%s: %s" % (cls, e.args[0])

        if log.getEffectiveLevel() > logging.DEBUG:     # no debug
            detail = production_error_msg % \
                        rnaseqlyze.admin_email
            body_template += "${detail}"
        else:                                           # debug
            detail = ''
            if isinstance(e, DBAPIError):
                detail += dberror_msg
            import traceback
            detail += '%s\n\nStack trace:\n' % e
            detail += ''.join(traceback.format_tb(exc_info[2]))

            log.debug(detail)
            body_template += "<pre>\n${detail}</pre>"

        HTTPError.__init__(self, detail, body_template=body_template)

dberror_msg = """\
This is a database related eror.

If it is not yet initialized or the schema has changed,
just run the "rnas-dbinit" script to (re-)initialize it.

Afterwards, restart the Pyramid application, i.e. send a
SIG_INT to the apache mod_wsgi daemon processes, and try again.
"""

production_error_msg = """\
If you think that this is a bug, please contact the application administrator,
%s, and inform him/her of the time the error occurred, and anything you might
have done that may have caused the error.
Thank You!
"""
# 'You' is intentionally capitalized! :-) Rule 84: http://goo.gl/BLBwX

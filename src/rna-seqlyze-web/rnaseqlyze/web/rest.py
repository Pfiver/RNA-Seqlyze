"""
Pyramid REST Views
"""

import logging
log = logging.getLogger(__name__)

from pyramid.view import view_config

import rnaseqlyze
from rnaseqlyze.web import DBSession, DBSession_unmanaged
from rnaseqlyze.core import service
from rnaseqlyze.core.orm import Analysis

@view_config(route_name='analysis_rest', renderer='json')
def display_rest(request):
    """
    **REST Analysis View**
    """
    id = int(request.matchdict["id"])
    items = vars(DBSession.query(Analysis).get(id)).items()
    return dict((k, str(v)) for (k, v) in items if k[0] is not '_')

@view_config(route_name='analysis_files_rest', renderer='json')
def analysis_files_rest(request):
    """
    **REST Files View**

    This view provides a (minimalistic, only GET is
    implemented) REST interface to '/analysis/{id}/files'.
    """
    return [ { "path": "foo" }, { "path": "bar" } ]

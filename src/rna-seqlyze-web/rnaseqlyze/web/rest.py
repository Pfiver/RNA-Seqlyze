"""
Pyramid REST Views
"""
import logging
log = logging.getLogger(__name__)

import os

from pyramid.view import view_config

import rnaseqlyze
from rnaseqlyze.web import DBSession, DBSession_unmanaged
from rnaseqlyze.core import service
from rnaseqlyze.core.orm import Analysis

none_type = type(None)
def json_filter_attributes(kv):
    if kv[0][0] != '_' and \
       type(kv[1]) in (none_type, bool, int, long, float, str, list):
        return True
    return False

@view_config(route_name='analysis_rest', renderer='json')
def display_rest(request):
    """
    **REST Analysis View**
    """
    analysis = DBSession.query(Analysis).get(int(request.matchdict["id"]))
    return dict(filter(json_filter_attributes, analysis.__dict__.iteritems()))

@view_config(route_name='analysis_files_rest', renderer='json')
def analysis_files_rest(request):
    """
    **REST Files View**

    This view provides a (minimalistic, only GET is
    implemented) REST interface to '/analysis/{id}/files'.
    """
    files = []
    analysis = DBSession.query(Analysis).get(int(request.matchdict["id"]))
    os.chdir(analysis.data_dir)
    for dirpath, dirnames, filenames in os.walk("."):
        dir = dirpath[2:]
        for fn in filenames:
            files.append({'path': os.path.join(dir, fn)})
    return files

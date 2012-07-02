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
from rnaseqlyze.core.orm import Analysis, UCSCOrganism

@view_config(route_name='analysis_rest', renderer='jsonx')
def display(request):
    """
    **REST Analysis View**
    """
    return DBSession.query(Analysis).get(int(request.matchdict["id"]))

@view_config(route_name='analysis_files_rest', renderer='jsonx')
def analysis_files(request):
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

@view_config(route_name='organisms_rest', renderer='jsonx')
def organisms(request):
    """
    ***REST Organisms View***

    Displays the list of organism titles
    along with their UCSC db and NCBI RefSeq accession identifiers
    """
    return DBSession.query(UCSCOrganism).all()

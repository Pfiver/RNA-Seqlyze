"""
Tools to deal with the UCSC genome browser at http://archaea.ucsc.edu/
"""
import logging
log = logging.getLogger(__name__)

from json import load
from urllib2 import urlopen
from urlparse import urljoin
from StringIO import StringIO
from shutil import copyfileobj
from os import listdir, makedirs
from os.path import join, dirname, isdir

from lxml.html import parse
from lxml.etree import dump

import rnaseqlyze
from rnaseqlyze.core import security
from rnaseqlyze.core.orm import UCSCOrganism

cart_reset_url = "http://archaea.ucsc.edu/cgi-bin/cartReset"
custom_track_url = "http://archaea.ucsc.edu/cgi-bin/hgTracks"
custom_track_params = "?db={org_db}&hgt.customText={track_url}"

# FIXME: This part need some work:
# 1) The module presently can only be imported after
#    rnaseqlyze.configure(workdir) was called.
# 2) The org_list_default_dir = dirname(__file__)
#    hack will not work if the distribution is installed
#    as a zipped .egg. pkg_resources.resource_stream or
#    pkg_resources.resource_string should be used instead.
org_list_base_url = "http://archaea.ucsc.edu/wp-content/data/"
org_list_cache_dir = join(rnaseqlyze.workdir, "ucsc-orglist-cache")
org_list_default_dir = dirname(__file__)
json_links_file_name = "ucsc-wp-data.html"


def get_org_list():

    if not isdir(org_list_cache_dir):
        makedirs(org_list_cache_dir)

    return list(get_organisms(get_json_files()))

def get_json_files():

    json_links_file = None
    json_files = None

    for get_json_links_file in (get_json_links_file_web,
                                get_json_links_file_cache,
                                get_json_links_file_default,):
        try:
            log.debug("trying %s" % get_json_links_file.func_name)
            json_links_file = get_json_links_file()
        except Exception, e:
            log.warn("%s failed: %r" % (get_json_links_file.func_name, e))
            continue
        
        for get_json_files in (get_get_json_files_web(json_links_file),
                               get_json_files_cache,
                               get_json_files_default,):

            try:
                log.debug("trying %s" % get_json_files.func_name)
                return get_json_files()
            except Exception, e:
                log.warn("%s failed: %r" % (get_json_files.func_name, e))
    
    raise Exception("Couldn't get json organism lists")

# getting the links file

def get_json_links_file_web():
    json_links_file = StringIO()
    copyfileobj(urlopen(org_list_base_url), json_links_file)
    json_links_file.seek(0)
    return json_links_file

def get_json_links_file_cache():
    return open(join(org_list_cache_dir, json_links_file_name))

def get_json_links_file_default():
    return open(join(org_list_default_dir, json_links_file_name))


# getting json files

def get_get_json_files_web(json_links_file):

    def get_json_files_web():
        for e in parse(json_links_file).getroot().iter("a"):
            link = e.attrib['href']
            if link.endswith(".json"):
                security.check_valid_filename(link)
                # FIXME: The json files should also be returned as StringIO
                #        buffers only and the cache files shouldn't be
                #        overwritten until it is certain that the newly
                #        downloaded files contain the expected data
                json_file = open(join(org_list_cache_dir, link), "w+")
                copyfileobj(urlopen(urljoin(
                        org_list_base_url, link)), json_file)
                json_file.seek(0)
                yield json_file

        try:
            json_links_file.fileno()
            # json_links_file was defaults or cached
        except:
            # json_links_file was a memory buffer - save it because
            # if this code is is reached, it means there was no error,
            # in the code above, so the buffer likely contains good links
            json_links_file.seek(0)
            copyfileobj(json_links_file, open(join(
                    org_list_cache_dir, json_links_file_name), "w"))

    return get_json_files_web

def get_json_files_cache():
    for json in listdir(org_list_cache_dir):
        if json.endswith(".json"):
            yield open(join(org_list_cache_dir, json))

def get_json_files_default():
    for json in listdir(org_list_default_dir):
        if json.endswith(".json"):
            yield open(join(org_list_default_dir, json))


# creating UCSCOrganism objects from json files

def get_organisms(json_files):
    for json_file in json_files:
        for object in load(json_file):
            for child in object['children']:
                for grandchild in child['children']:
                    if grandchild['attr']['rel'] == 'genome':
                        yield UCSCOrganism(db=grandchild['attr']['id'],
                                           title=grandchild['data']['title'])

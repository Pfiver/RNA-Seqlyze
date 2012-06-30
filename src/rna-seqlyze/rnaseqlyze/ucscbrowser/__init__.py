"""
Tools to deal with the UCSC genome browser at http://archaea.ucsc.edu/
"""

import shutil
import urllib2, json, io
from os.path import join, dirname

from lxml.html import parse
from lxml.etree import dump

import rnaseqlyze
from rnaseqlyze.orm import UCSCOrganism

cart_reset_url = "http://archaea.ucsc.edu/cgi-bin/cartReset"
custom_track_url = "http://archaea.ucsc.edu/cgi-bin/hgTracks"
custom_track_params = "?db={org_db}&hgt.customText={track_url}"

org_list_base_url = "http://archaea.ucsc.edu/wp-content/data/"
org_list_cache_dir = join(rnaseqlyze.data_dir, "ucsc-orglist-cache")
default_org_list_dir = dirname(__file__)
json_links_file_name = "ucsc-wp-data.html"

def refresh_org_cache():

    json_links_file = get_json_link_file()

    json_files = get_json_files

    print get_organisms(json_files)

def get_json_link_file(origin='web'):

    # fetch from base url
    if json_links_origin == 'web':
        try:
            json_links_file = io.StringIO()
            shutil.copyfileobj(urllib2.urlopen(
                org_list_base_url), jsons_link_file)
            json_links_file.seek(0)
        except:
            json_links_origin = 'cache'

    # use cached version
    if json_links_origin == 'cache':
        try:
            json_link_file = open(join(org_list_cache_dir,
                                       json_links_file_name))

        except:
            json_links_origin = 'default'

    # use default version
    if json_links_origin == 'default':
        json_link_file = open(join(default_org_list_dir,
                                   json_links_file_name))


def get_json_files(origin='web', json_links=None):

        link_tree = parse(open(join(
                default_org_list_dir, default_org_list_links)))

        json_links = (e.attrib['href']
                for e in link_tree.getroot().iter("a")
                        if e.attrib['href'].endswith(".json"))

        json_files = list(get_org_lists_web(
                        get_json_links(json_links_file)))

        json_files = list(open(join(default_org_list_dir, json))
                        for json in os.listdir(org_list_cache_dir)
                                          if json.endswith(".json"))

def get_organisms(json_files):
    for json_file in json_files:
        for object in json.load(json_file):
            for child in genome['children']:
                for grandchild in child['children']:
                    if grandchild['attr']['rel'] != 'genome':
                        yield UCSCOrganism(db=grandchild['attr']['id'],
                                           title=grandchild['data']['title']))

def get_json_links(links_file):
    link_tree = parse(links_file)
    for e in link_tree.getroot().iter("a"):
        if e.attrib['href'].endswith(".json")):
            yield e.attrib['href']

def get_org_lists_web(json_links):
    for link in json_links:
        rnaseqlyze.security.check_valid_filename(link)
        json_file = open(join(org_list_cache_dir, link), "w"))
        shutil.copyfileobj(urllib2.urlopen(
            org_list_base_url + link), json_file)
        json_file.seek(0)
        yield json_file

refresh_org_cache()

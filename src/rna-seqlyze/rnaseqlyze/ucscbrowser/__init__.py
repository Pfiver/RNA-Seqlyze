"""
Tools to deal with the UCSC genome browser at http://archaea.ucsc.edu/
"""

import urllib2, json
from os.path import join, dirname

from lxml.html import parse
from lxml.etree import dump

import rnaseqlyze
from rnaseqlyze.orm import UCSCOrganism

cart_reset_url = "http://archaea.ucsc.edu/cgi-bin/cartReset"
custom_track_url = "http://archaea.ucsc.edu/cgi-bin/hgTracks" \
                            + "?db={db}&hgt.customText={customText}"

org_list_base_url = "http://archaea.ucsc.edu/wp-content/data/"

org_list_cache_dir = join(rnaseqlyze.data_dir, "ucsc-orglist-cache")

default_org_list_dir = dirname(__file__)
default_org_list_links = "ucsc-wp-data.html"

def refresh_org_cache():

    error = None
    organisms = []
    json_files = []

    try:
        tree = parse(urllib2.urlopen(org_list_base_url))
    except Exception, e:
        log.error("Couldn't retrieve organism list links.")
        log.error("The error was: %r" % e)
        log.info("Using defaults in %s" % default_org_list_dir)
        tree = parse(open(join(default_org_list_dir, default_org_list_links)))

    json_links = (e.attrib['href']
            for e in tree.getroot().iter("a")
                    if e.attrib['href'].endswith(".json"))

        try:
            for link in json_links:
                rnaseqlyze.security.check_valid_filename(link)
                json_file = open(join(org_list_cache_dir, link), "w"))
                shutil.copyfileobj(urllib2.urlopen(
                    org_list_base_url + link), json_file)
                json_files.append(json_file)
                json_file.seek(0)
        except Exception, e:
            error = e
            log.error("Couldn't retrieve json organism lists.")
            log.error("The error was: %r" % e)
            log.info("Using defaults in %s" % default_org_list_dir)
            for json in os.listdir(org_list_cache_dir):
                if json.endswith(".json"):
                    json_files.append(open(join(default_org_list_dir, json)))
    else:
        for link in json_links:
            json_files.append(open(join(default_org_list_dir, link)))

    for json_file in json_files:
        for object in json.load(json_file):
            for child in genome['children']:
                for grandchild in child['children']:
                    if grandchild['attr']['rel'] != 'genome':
                        continue
                    organisms.append(
                            UCSCOrganism(db=grandchild['attr']['id'],
                                         title=grandchild['data']['title']))

    print organisms
    
refresh_org_cache()

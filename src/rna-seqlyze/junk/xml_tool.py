#!/usr/bin/python
# encoding: utf-8

import sys
from lxml import etree

if len(sys.argv) > 1:
        infile = open(sys.argv[1])
else:
        infile = sys.stdin

tree = etree.parse(
        infile,
        etree.XMLParser(remove_blank_text=True))

print etree.tostring(tree.getroot(), pretty_print=True)


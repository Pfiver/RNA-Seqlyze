"""\
RNA-Seqlyze xmltool

XML version of of `python -m json.tool`.

Takes an xml file or stream as input and pretty-prints it.

Usage:
    rnas-xmltool <input.xml> <output.xml>

If <input.xml> is '-', use 'sys.stdin, if <output.xml> is '-', use 'sys.stdout'.
"""
import sys
from lxml import etree

def main():

    if len(sys.argv) < 2 or sys.argv[1] in ('-h', '--help'):
        print __doc__
        return

    inputfile = sys.argv[1] == '-' and sys.stdin or open(sys.argv[1])
    outputfile = sys.argv[2] == '-' and sys.stdout or open(sys.argv[2], "w")

    tree = etree.parse(inputfile,
                       etree.XMLParser(remove_blank_text=True))
    print >> outputfile, etree.tostring(tree.getroot(), pretty_print=True)

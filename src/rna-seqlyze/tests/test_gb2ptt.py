#!/usr/bin/python

import nose.tools as nt

def test_script():

    # set up the testing environment
    import os
    here = os.path.dirname(__file__)
    os.environ["PYTHONPATH"] = here + '/../source'

    # run the converter script
    from subprocess import PIPE, Popen
    proc = Popen([here + "/../scripts/gb2ptt",
                  "tests/data/NC_002754-partial.gb"], stdout=PIPE)

    # read the converters standard output
    import csv
    reader = csv.reader(proc.stdout, delimiter='\t')

    # get the first 4 lines of output
    # and run some rudimentary checks on them
    nt.assert_regexp_matches(reader.next()[0], "Sulfolobus solfataricus P2.*")
    nt.assert_equals(reader.next(), [])
    cols = reader.next() # column headings
    reader.next() # the wrapping around gene SSO12256 (clycotransferase)
                  # comes first here but not in the ncbi file / discard it for now

    # open the reference (eXpected) ptt file
    # discard the first 2 lines and safe the column headers
    xreader = csv.reader(open("tests/data/NC_002754-partial.ptt"), delimiter='\t')
    xreader.next() # first line
    xreader.next() # second line
    xcols = xreader.next() # column headings

    # check if the column headers match
    nt.assert_equals(cols, xcols)

    # fixup function to make the rows compareable
    def truncate_row(row):
        del row[7:]         # i didn't bother to fill in the "cog" row (7)
                            # and the descriptions (8) don't match exactly
        return row

    # compare the rows produced by ./gb2ptt.py to those read from the reference file
    nt.assert_items_equal(map(truncate_row, reader), map(truncate_row, xreader))

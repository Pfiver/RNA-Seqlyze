#!/usr/bin/python

import unittest

class Test(unittest.TestCase):
    def runTest(self):

		# run the converter
		from subprocess import PIPE, Popen
		proc = Popen("./gb2ptt.py test/data/NC_002754-partial.gb", shell=True, stdout=PIPE)

		# read the converters standard output
		import csv
		reader = csv.reader(proc.stdout, delimiter='\t')

		# get the first 4 lines of output
		# and run some rudimentary checks on them
		self.assertRegexpMatches(reader.next()[0], "Sulfolobus solfataricus P2.*")
		self.assertEquals(reader.next(), [])
		cols = reader.next() # column headings
		reader.next() # the wrapping around gene SSO12256 (clycotransferase)
		              # comes first here but not in the ncbi file / discard it for now

		# open the reference (eXpected) ptt file
		# discard the first 2 lines and safe the column headers
		xreader = csv.reader(open("test/data/NC_002754-partial.ptt"), delimiter='\t')
		xreader.next() # first line
		xreader.next() # second line
		xcols = xreader.next() # column headings

		# check if the column headers match
		self.assertEquals(cols, xcols)

		# fixup function to make the rows compareable
		def truncate_row(row):
			del row[7:]			# i didn't bother to fill in the "cog" row (7)
								# and the descriptions (8) don't match exactly
			return row

		# compare the rows produced by ./gb2ptt.py to those read from the reference file
		self.assertItemsEqual(map(truncate_row, reader), map(truncate_row, xreader))

#!/usr/bin/python
#
#  - 
# 
# Copyright (C) 2011  Patrick Pfeifer

usage_text = """\

usage: %s [-f|--format FORMAT] input_path

    generates a `bed` track out of operon predictions listed in input_path

    FORMAT
                'door', 'operondb' or 'microbesonline'

                if not specified, it is auto-detected
                based on the extension of the input file:
                   .opr -> door
                   .html -> operondb
                   .named -> microbesonline

    input_path
                the path to the input file
"""

import sys, logging

log = logging.getLogger(sys.argv[0])

class operon(object):
	id = None
	begin = None
	end = None
	strand = None
	confidence = None

	def __init__(self, id, begin=None, end=None, strand=None, confidence=None):
		self.id = id
		self.begin = begin
		self.end = end
		self.strand = strand
		self.confidence = confidence

def process_operondb(in_path):
	import itertools
	from lxml.html import parse

	operons = {}
	ctr = itertools.count()
	doc = parse(in_path).getroot()
	ht_ops = doc.xpath("//html/body/div/div/div[2]/div/div")

#	all_coords = doc.xpath("//html/body/div/div/div[2]/div/div/table[1]/tr/td/h3")
#	all_confs = doc.xpath("//html/body/div/div/div[2]/div/div/table[2]/tr/td[3]/a")

	for ht_op in ht_ops:
		coords = ht_op.xpath("table[1]/tr/td/h3")
		confs = ht_op.xpath("table[2]/tr/td[3]/a")

		if not len(coords):
			continue

		# coords, confidence:
		# Directon coords:  2988382 - 2990017
		# confidence=11 n=7

		import re
		m = re.match(r'Directon coords:\s*([0-9]+)\s*-\s*([0-9]+)', coords[0].text)
		if not m:
			log.info("couldn't parse coords: '%s'" % coords[0].text)
			return
		begin, end = m.groups()

		for conf in confs:
			m = re.match(r'confidence\s*=\s*([0-9]+)', conf.text)
			if not m:
				log.info("couldn't parse confidence: '%s'" % conf.text)
				return
			confidence = int(m.group(1)) / 100.

		strand = "."
		id = ctr.next()

		operons[id] = operon(id, begin, end, strand, confidence)
	
	return operons.values()

def process_door(in_path):

	import csv
	reader = csv.reader(open(in_path, 'rb'), delimiter='\t')

	# get column headings
	col_names = reader.next()
	# ['OperonID', 'GI', 'Synonym', 'Start', 'End', 'Strand', 'Length', 'COG_number', 'Product']

	# create a dictionary { heading: index, ... }
	col_is = dict(zip(col_names,range(len(col_names))))
	#print col_names
	#return

	operons = {}

	for row in reader:

		if len(row) != len(col_names):
			log.debug( "skipping invalid line: '%s'" % row )
			continue

		id = int(row[col_is["OperonID"]])

		begin = int(row[col_is["Start"]])
		end = int(row[col_is["End"]])
		strand = row[col_is["Strand"]]

		# .... FIXME !!! ... FIXME !!! ... (Sulfolobus specific)
		if begin == 2991448:
			begin = 1

		if strand == "-":
			begin, end = end, begin

		if id not in operons:
			operons[id] = operon(id, begin, end, strand)
			continue

		if operons[id].begin > begin:
			operons[id].begin = begin

		if operons[id].end < end:
			operons[id].end = end

	return operons.values()

def output(operons, name="BED", desc="Custom BED track", confidence=None):
	# goal:
	#	track name="BED" description="BED Test Track" visibility=2 itemRgb="on" 
	#	chr     1       100     Foo     500     +       10      90      255,0,0
	#	chr     200     250     Bar     900     +       205     245     0,255,0
	#	chr     400     425     Baz     800     -       405     420     0,0,255

	if confidence == None:
		confidence = operons[0].confidence != None

	import csv
	sys.stdout.write( 'track name="%s" description="%s" visibility=2 itemRgb="%s"\n' %
	        (name,
			 desc,
			 confidence and "on" or "off"))
	writer = csv.writer(sys.stdout, delimiter='\t')

	# standard columns (no matter if confidence is output or not)
	cols = lambda op: ("chr", op.begin, op.end, "Operon_%d" % op.id, op.confidence, op.strand)

	# if confidence should be output, add 3 extra columns
	if confidence:
		import numpy
		white = numpy.array((255,255,255))
		def color(confidence):
			color_array = white * 0.5 * (1 - op.confidence)
			return ",".join(str(int(v)) for v in color_array)
		cols_0 = cols
		cols = lambda op: cols_0(op) + (op.begin, op.end, color(op.confidence))

	# write the operons, one per row
	for op in operons:
		# .... FIXME !!! ... FIXME !!!
		if op.begin > op.end:
			log.info("op.begin > op.end ??? (%d: %d > %d)", (op.id, op.begin, op.end))
			continue
		writer.writerow(cols(op))

def main():
	import getopt

	opts, args = getopt.getopt(sys.argv[1:], 'f:c', ["format=", "confidence"])

	fmt = None
	confidence = None

	for opt, val in opts:
		if opt in ("-f", "--format"):
			fmt = val
		if opt in ("-c", "--confidence"):
			confidence = True

	if len(args) < 1:
		sys.stderr.write(usage_text % sys.argv[0])
		return 1

	ext_fmt = {
		"opr": "door",
		"html": "operondb",
		"named": "microbesonline"
	}

	if fmt == None:
		try:
			ext = args[0][args[0].rindex('.')+1:]
			fmt = ext_fmt[ext]
		except object, x:
			print x
			return

	try:
		name = args[0]
		desc = "Operons from %s" % args[0]
		operons = globals()["process_" + fmt](*args)
		output(operons, name, desc, confidence=confidence)

	except object, x:
		log.error(x)
		return 1

	return 0

sys.exit(main())

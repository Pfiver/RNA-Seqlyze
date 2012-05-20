#!/usr/bin/env python
# encoding: utf-8

# this "distrubute" setup script depends on

	# "distribute" 0.6.26, by Tarek Ziadé
	# - source code: https://bitbucket.org/tarek/distribute
	# - documentation: http://pypi.python.org/pypi/distribute

	# to run this script, "distribute" must be installed, as follows:
	# 1) fetch the setup script
	#     $ curl -O https://bitbucket.org/tarek/distribute/raw/b69f072c0002/distribute_setup.py
	# 2) install "distribute"
	#     $ python distribute_setup.py --user

from setuptools import setup
from setuptools import find_packages

setup(
	name='RNA-seqlyze',
	version='0.1',
	test_suite='test',
	license = "Mixed",
	author_email = "patrick@patrickpfeifer.net"
	description = "RNA-seq analysis & sequence annotation enhancement web-application",
	license = "MIT",
	url = "http://git.pfeifer.ch/rna-seqlyze.git",
	packages = find_packages(),
    scripts = [
		'scripts/xml_tool',
	],
	install_requires = [
#		"Cython >= 0.16",
		"pysam >= 0.6",
		"Trac >= 0.13dev",
		"MySQL-python 1.2.2",	# for Trac wiki
		"MarkupSafe >= 0.15",	# required by Trac
		"biopython >= 1.59",
		"bcbio_nextgen >= 0.3a",
	],
)

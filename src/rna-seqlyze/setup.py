#!/usr/bin/python
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

setup(
	name='RNA-seqlyze',
	version='0.1',
	test_suite='test'
)

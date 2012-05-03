#!/usr/bin/python

import distribute_setup
distribute_setup.use_setuptools()
from setuptools import setup

setup(
	name='RNA-seqlyze',
	version='0.1',
	test_suite='test')

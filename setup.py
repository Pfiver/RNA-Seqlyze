#!/usr/bin/python

from distutils.core import Command, setup
from unittest import TextTestRunner, defaultTestLoader

# custom setup commands:
#  https://da44en.wordpress.com/2002/11/22/using-distutils/

class TestCommand(Command):
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        TextTestRunner().run(
			defaultTestLoader.discover("test", "*"))

setup(
	name='RNA-seqlyze',
	version='0.1',
	cmdclass={
		'test': TestCommand
	})

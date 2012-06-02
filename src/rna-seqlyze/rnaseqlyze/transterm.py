#!/usr/bin/python

def run(*args):

	from os.path import abspath, dirname
	tt_dir = dirname(abspath(__file__)) + "/../../TransTermHP/"
	tt_cmd = tt_dir + "transterm", "-p", tt_dir + "expterm.dat"

	import subprocess
	subprocess.call(tt_cmd + args)

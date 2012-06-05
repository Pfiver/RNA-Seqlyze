# encoding: utf-8

def run(*args):

	from os.path import abspath, dirname
	tt_dir = dirname(abspath(__file__)) + "/../../TransTermHP/"
	tt_cmd = tt_dir + "transterm", "-p", tt_dir + "expterm.dat"

	import subprocess
    if subprocess.call(tt_cmd + args) != 0:
        raise Exception(str(tt_cmd + args) + " failed")

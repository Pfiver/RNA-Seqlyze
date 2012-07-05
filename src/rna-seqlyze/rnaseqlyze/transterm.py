"""
A module to run transterm_hp
"""

import os, subprocess
from subprocess import PIPE

def run(args, out=None, err=None):
    """
    Run transterm_hp with the given arguments plus "-p expterm.dat"
    """
    def findit():
        for path in os.getenv("PATH").split(os.path.pathsep):
            for name in os.listdir(os.path.join(path, "../lib")):
                if name == 'expterm.dat':
                    return os.path.join(path, "../lib", name)
    expterm_dat = findit()
    if not expterm_dat:
        raise Exception("'expterm.dat' not found")
    cmd = 'transterm', '-p', expterm_dat
    proc = subprocess.Popen(cmd + args, stdout=out, stderr=err)
    proc.wait()
    if proc.returncode != 0:
        raise Exception(str(cmd + args) + " failed")

def fa2bed(tt_output, bed_file):

    for line in tt_output:

        if not line.startswith("  TERM"):
            return

        term, id, beg, dash, end, str, pos, con, rest = line.split(None, 9)

        col = 100 - int(con)    # let color vary from 0 (black) to 100 (gray)
        if str == '-':          # switch begin & end if on reverse strand
            beg, end = end, beg

        print >> bed_file, '\t'.join(('chr', beg, end,
                                      'TERM_' + id, con, str,
                                      beg, end, ','.join((str(col),)*3)))

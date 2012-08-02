"""
A module to run transterm_hp
"""

import os, subprocess

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
    cmd = ('transterm', '-p', expterm_dat) + tuple(args)
    proc = subprocess.Popen(cmd, stdout=out, stderr=err)
    proc.wait()
    if proc.returncode != 0:
        raise Exception(str(cmd) + " failed")

def tt2bed(tt_output, bed_file):
    for id, begin, end, strand, confidence in iterator(tt_output):
        # let color vary from 0 (black) to 100 (gray)
        rgb_color = ','.join((str(100 - int(confidence)),)*3)
        print >> bed_file, '\t'.join((
            'chr', begin, end, 'TERM_' + id,
            str(confidence), strand, begin, end, rgb_color
        ))

def iterator(tt_output):
    for line in tt_output:
        if not line.startswith("  TERM"):
            continue
        TERM, id, begin, dash, end, \
              strand, position, confidence, rest = line.split(None, 8)
        # switch begin & end on reverse strand
        if strand == '-':
            begin, end = end, begin
        yield id, begin, end, strand, int(confidence)

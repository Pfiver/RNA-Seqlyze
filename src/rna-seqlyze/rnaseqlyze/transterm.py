"""
run transterm_hp with the given arguments plus "-p expterm.dat"
"""


def run(args):
    import os, subprocess
    def findit():
        for path in os.getenv("PATH").split(os.path.pathsep):
            for name in os.listdir(os.path.join(path, "../lib")):
                if name == 'expterm.dat':
                    return os.path.join(path, "../lib", name)
    expterm_dat = findit()
    if not expterm_dat:
        raise Exception("'expterm.dat' not found")
    cmd = ['transterm', '-p', expterm_dat]
    proc = subprocess.Popen(cmd + args)
    if proc.returncode != 0:
        raise Exception(str(cmd + args) + " failed")

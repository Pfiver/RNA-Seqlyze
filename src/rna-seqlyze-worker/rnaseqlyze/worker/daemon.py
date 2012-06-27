"""
RNA-Seqlyze Worker

Usage:
    rnas-worker (start|stop|restart)
    rnas-worker --development

Start, stop, restart the worker daemon.
"""

def main():
    import os
    import docopt
    import paste.script.exe

    opts = docopt.docopt(__doc__)
    for command in "start|stop|restart".split('|'):
        if opts[command]:
            mode = "production"
            args = [command]
            break
    else:
        mode = "development"
        args = []

    if mode == 'production':
        import rnaseqlyze
        args.extend([
            "--user=" + rnaseqlyze.worker_user,
            "--group=" + rnaseqlyze.group,
        ])

    here = os.path.dirname(os.path.abspath(__file__))
    conf_file = os.path.join(here, '..', '..', mode + '.ini')

    args.insert(0, conf_file)
    os.environ['_'] = conf_file
    paste.script.exe.ExeCommand("exe").run(args)

"""
RNA-Seqlyze Worker

Usage:
    rnas-worker (start|stop|restart)
    rnas-worker --development

Start, stop or restart the worker daemon.
"""

def main():
    import os
    import docopt
    import paste.script.serve

    opts = docopt.docopt(__doc__)
    for command in "start|stop|restart".split('|'):
        if opts[command]:
            mode = "production"
            args = [command, "--daemon"]
            break
    else:
        mode = "development"
        args = ["--reload"]

    if mode == 'production':
        import rnaseqlyze
        args.extend([
            "--user=" + rnaseqlyze.worker_user,
            "--group=" + rnaseqlyze.group,
            "--log-file=" + rnaseqlyze.worker_log_file,
            "--pid-file=" + rnaseqlyze.worker_pid_file,
        ])

    here = os.path.dirname(os.path.abspath(__file__))
    conf_file = os.path.join(here, '..', '..', mode + '.ini')

    args.insert(0, conf_file)
    paste.script.serve.ServeCommand("serve").run(args)

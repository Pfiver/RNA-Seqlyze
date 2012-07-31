"""
RNA-Seqlyze Worker

Start, stop or restart the worker daemon
or run it in the foreground, in development mode.

Usage:
    rnas-worker <workdir> (start|stop|restart)
    rnas-worker <workdir> --development
    rnas-worker -h|--help

Arguments:

    <workdir> The path to the workers 'workdir'.
              The 'workdir' is where the configuration, the
              application database and all analysis data are stored.

    start|stop|restart

              If one of those arguments is given, the daemon is
              run in the background. It will write it's PID to the
              file <workdir>/worker-daemon.pid and its output will be logged
              to <workdir>/worker-daemon.log. This is not the "log file"
              however. The "log file" path can be configured
              in <workdir>/worker.ini.

    --development

              If this argument is present, the worker daemon is run in
              development mode, which means that it will no fork to the
              background. If any source files (.py) are changed when the
              daemon is running in development mode, it will be
              automatically restarted.

"""

from os.path import abspath, join

from paste.script import serve
import docopt

import rnaseqlyze

def main():
    opts = docopt.docopt(__doc__)

    for command in "start|stop|restart".split('|'):
        if opts[command]:
            mode = "production"
            args = [command, "--daemon"]
            break
    else:
        mode = "development"
        args = ["--reload"]

    workdir = abspath(opts['<workdir>'])
    rnaseqlyze.configure(workdir)

    if mode == 'production':
        args.extend([
            "--log-file=" + join(workdir, 'worker-daemon.log'),
            "--pid-file=" + join(workdir, 'worker-daemon.pid'),
        ])
    args.extend([
        "worker_port=" + rnaseqlyze.worker_port,
    ])

    conf_file = join(workdir, 'worker.ini')
    serve.ServeCommand("serve").run([conf_file] + args)


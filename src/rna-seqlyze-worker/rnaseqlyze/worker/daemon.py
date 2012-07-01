"""
RNA-Seqlyze Worker

Start, stop or restart the worker daemon
or run it in the foreground, in development mode.

Usage:
    rnas-worker <workdir> (start|stop|restart)
    rnas-worker <workdir> --development
    rnas-worker -h|--help

Arguments:
    <workdir>     The path to the workers 'workdir'.
                  The 'workdir' is where the configuration, the
                  application database and all analysis data are stored.

    start         If one of those arguments is given, the daemon is
    stop          run in the background. It will write it's PID to the
    restart       file <workdir>/worker-daemon.pid and its output will be logged
                  to <workdir>/worker-daemon.log. This is not the "log file"
                  however. The "log file" path can be configured
                  in <workdir>/worker.ini.

    --development If this argument is present, the worker daemon is run in
                  development mode, which means that it will no fork to the
                  background. If any source files (.py) are changed when the
                  daemon is running in development mode, it will be
                  automatically restarted.

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
            "--log-file=" + here + os.sep + 'worker-daemon.log',
            "--pid-file=" + here + os.sep + 'worker-daemon.pid',
        ])

    conf_file = os.path.join(opts['<workdir>'], 'worker.ini')
    paste.script.serve.ServeCommand("serve").run([conf_file] + args)


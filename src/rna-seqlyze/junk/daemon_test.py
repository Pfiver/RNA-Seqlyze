import daemon
with daemon.DaemonContext():
    import time
    import logging
    logging.basicConfig(filename="/home/pfeifer/daemon_test.log", level=logging.NOTSET)
    log = logging.getLogger("daemon_test")
    import os, pprint
    log.info(pprint.pformat(os.environ))
    for i in range(10):
        log.info("It's " + time.strftime("%Y-%m-%dT%H:%M:%S"))
        time.sleep(1)

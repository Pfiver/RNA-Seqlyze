#!/bin/bash

### BEGIN INIT INFO
# Provides:          rna-seqlyze-worker
# Required-Start:    $local_fs $network
# Required-Stop:     $local_fs $network
# Default-Start:     2 3 4 5
# Default-Stop:      0 1 6
# Short-Description: RNA-Seqlyze worker daemon
# Description:       Debian init script for the RNA-Seqlyze worker daemon
### END INIT INFO

# Usage:
#	- configure the two directories below
#       - copy this file to /etc/init.d/rnaseqlyze.sh
#	- enable it with `insserv /etc/init.d/rnaseqlyze.sh`

###

PREFIX=/home/biopython/.local		# the <prefix> passed to `rnas-build`
WORKDIR=/home/biopython/rnas-workdir	# the <workdir> passed to `rnas-init`

###

test -x $PREFIX/bin/rnas-worker || exit 0

PATH=$PREFIX/bin:/bin:/usr/bin

. /lib/lsb/init-functions

case "$1" in
    start|stop|restart)
        log_daemon_msg "${1^} RNA-Seqlyze worker daemon" rna-seqlyze-worker
	rnas-worker $WORKDIR $1
        log_end_msg $?
    ;;
    *)
        echo "Usage: /etc/init.d/rnaseqlyze.sh {start|stop|restart}"
        exit 1
    ;;
esac

exit 0

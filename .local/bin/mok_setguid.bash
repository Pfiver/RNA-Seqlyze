#!/bin/bash -e
exec 77<>$(mktemp)
fd=/proc/self/fd/77
rm -f $(readlink $fd)
echo 'setgid(){return 0;}int(*setuid)()=setgid;' |
        gcc -shared -o $fd -x c -; LD_PRELOAD=$fd exec "$@"

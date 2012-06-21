from __future__ import print_function

import sys
def main(argv=sys.argv):
    import pyinotify
    mask = sum(getattr(pyinotify, "IN_" + m) for m in
               "CREATE DELETE MOVED_TO MODIFY".split())
    from pyinotify import WatchManager, ProcessEvent, Notifier
    wm = WatchManager()
    for d in argv[1:]:
        wm.add_watch(d, mask, rec=True)
    class PrintEvents(ProcessEvent):
        def process_default(self, event):
            print(event.maskname, event.name)
            sys.stdout.flush()
    eh = PrintEvents()
    Notifier(wm, eh).loop()

"""\
Programs runnable from the command line.

For each module contained in this package, a wrapper script
called `rnas-<module name>` will be installed in `<prefix>/bin`.
"""

from .. import project_name
project_name += "-cli"

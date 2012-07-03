"""
Top level package module

Importing this package configures the python "logging"
module in a way that messages of any level go to sys.stderr.
"""

#: the project base-name
#: the -cli, -web and -worker project names
#: are constructed by appending the part name to this one
project_name = "rna-seqlyze"

import pkg_resources
#: The __version__ property is set automatically set to the value of
#: pkg_resources.get_distribution(project_name).version on module import time.
try:
    __version__ = pkg_resources.get_distribution(project_name).version
except:
    # When the packe is initially installed, this module is imported by setup.py
    # and the the `project_name` attribute defined above is used to set the
    # project name, in order to make that information non-redundant.
    #
    # The version, however, is determined at install/build time by running
    # `git --describe` - in setup.py. So it doesn't matter if it is not set at
    # that time.  Later on, because `python setup.py install/develop` writes it
    # to rna-seqlyze.egg-info/PKG-INFO, where pkg_resources picks it up,
    # it will be available whenever the package is imported.
    pass
del pkg_resources

import logging
logging.basicConfig(level=0, format="%(levelname)-5.5s [%(name)s] %(message)s")
del logging

def configure(_workdir):
    """
    Calling this function

        - sets rnaseqlyze.workdir to <workdir>
        
        - sets rnaseqlyze.<setting> attributes for all
          settings under [rnaseqlyze] in '<workdir>/rnaseqlyze.ini'.

        - imports Bio.Entrez and sets Bio.Entrez.email to rnaseqlyze.admin_email
    """

    global workdir
    workdir = _workdir

    from os.path import join
    from ConfigParser import ConfigParser
    config = ConfigParser(dict(here=workdir))
    config.read(join(workdir, 'rnaseqlyze.ini'))

    for name, value in config.items("rnaseqlyze"):
        globals()[name] = value

    import Bio.Entrez
    Bio.Entrez.email = admin_email

"""
Top level package module
"""

#: the project base-name
#: the -cli, -web and -worker project names
#: are constructed by appending the part name to this one
project_name = "rna-seqlyze"

import pkg_resources
#: The __version__ property is set automatically set to the value of
#: pkg_resources.get_distribution(project_name).version on module import time.
__version__ = pkg_resources.get_distribution(project_name).version
del pkg_resources

def configure(_workdir):
    """
    Calling this function

        - sets rnaseqlyze.workdir to <workdir>
        
        - sets rnaseqlyze.<setting> attributes for all
          settings under [rnaseqlyze] in '<workdir>/rnaseqlyze.ini'.

        - imports Bio.Entrez and sets Bio.Entrez.email to rnaseqlyze.admin_email

        - configures the python "logging" module by calling
          logging.config.fileConfig('<workdir>/rnaseqlyze.ini').
    """

    global workdir
    workdir = _workdir

    import os
    conf_ini = os.path.join(workdir, 'rnaseqlyze.ini')

    from ConfigParser import ConfigParser
    config = ConfigParser()
    config.read(conf_ini)

    for name, value in config.items("rnaseqlyze"):
        globals()[name] = value

    import Bio.Entrez
    Bio.Entrez.email = admin_email

    import logging.config.fileConfig
    logging.config.fileConfig(conf_ini)

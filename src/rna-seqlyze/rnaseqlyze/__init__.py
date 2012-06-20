"""
Top level package module
"""

project_name = 'rna-seqlyze'
config_filename = 'rnaseqlyze.ini'

def _init():
    import os.path
    import pkg_resources
    from ConfigParser import ConfigParser

    dist = pkg_resources.get_distribution(project_name)

    global __version__
    __version__ = dist.version

    global install_path
    install_path = dist.location

    config_path = os.path.join(install_path, config_filename)
    config = ConfigParser({'here': install_path})
    config.read(config_path)

    for name, value in config.items("rnaseqlyze"):
        globals()[name] = value

_init()
del _init

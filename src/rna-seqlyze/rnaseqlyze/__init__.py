project_name = 'rna-seqlyze'
config_filename = 'rnaseqlyze.ini'

def _get_conf():
    import os.path
    import ConfigParser
    import pkg_resources

    global install_path, config_path, config

    install_path = pkg_resources.resource_filename(
                   pkg_resources.Requirement.parse(project_name), '')

    config_path = os.path.join(install_path, config_filename)

    config = ConfigParser.ConfigParser({'here': install_path})
    config.read(config_path)

_get_conf()
del _get_conf

db_url = config.get("database", "url")
cache_path = config.get("cache", "path")

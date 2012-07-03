"""
RNA-Seqlyze WSGI Apllication

Provides the get_app(workdir) function,
which returns a wsgi application callable.
"""

def get_app(workdir):
    """
    Basically returns wrapper around paster.get_app
    that strips the ".wsgi" extension from SCRIPT_NAME
    """

    # configure the core package
    import rnaseqlyze
    rnaseqlyze.configure(workdir)

    # default configuraion file name
    import os.path
    web_ini = os.path.join(workdir, 'web.ini')

    # configure logging
    import logging.config
    logging.config.fileConfig(web_ini, dict(here=workdir))

    # create the pyramid wsgi app
    import pyramid.paster
    pyramid_app = pyramid.paster.get_app(web_ini, 'main')

    # return a wrapper that adjusts SCRIPT_NAME
    def app(environ, start_request):
        environ['SCRIPT_NAME'] = \
                environ['SCRIPT_NAME'][:-5]
        return pyramid_app(environ, start_request)

    return app

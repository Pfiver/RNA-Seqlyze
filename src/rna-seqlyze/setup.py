# rna-seqlyze setup script

    # depends on git and "distribute"
    # - git: git-scm.com
    # - "distribute":
    #  - docs: http://packages.python.org/distribute/
    #  - code: https://bitbucket.org/tarek/distribute
    #  - installation:
    #    $ curl -O https://bitbucket.org/tarek/distribute/raw/default/distribute_setup.py
    #    $ python distribute_setup.py --user

from setuptools import setup, find_packages

def get_version():
    from subprocess import Popen, PIPE
    git_describe = Popen(("git", "describe"), stdout=PIPE, stderr=PIPE)
    version, error = git_describe.communicate()
    if git_describe.returncode:
        raise Exception("couldn't determine package version: " + error)
    return version.strip()[1:]

setup(
    name='rna-seqlyze',
    version=get_version(),
    author="Patrick Pfeifer",
    author_email="patrick@patrickpfeifer.net",
    description="RNA-seq analysis | core features",
    long_description="RNA-seq analysis & sequence annotation enhancement web-application",
    license="Mixed",
    packages=find_packages(),
    test_suite='nose.collector',
    entry_points={
        'console_scripts': [
            'build-rnaseqlyze = rnaseqlyze.build:main',
            'rnaseqlyze-dbinit = rnaseqlyze.dbinit:main',
        ],
    },
    setup_requires=[
        "nose >= 1.0",
        "distribute >= 0.6.14",     # unsure if really required, but won't harm
    ],
    install_requires=[
        "pyflakes",                 # Buildbot
        "nose >= 1.1.2",
        "pyramid >= 1.3.2",
        "MarkupSafe >= 0.15",       # in Trac
        "MySQL-python >= 1.2.2",    # in Trac
    ],
)

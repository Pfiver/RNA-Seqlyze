from setuptools import setup, find_packages

import os, subprocess as sp
git_proc = sp.Popen("git describe --tags".split(), stdout=sp.PIPE)
git_out = git_proc.communicate()[0]
if git_proc.returncode:
    raise Exception()
version = git_out.strip()[1:]

setup(
    name=os.getcwd().split(os.sep)[-1],
    version=version,
    author="Patrick Pfeifer",
    author_email="patrick@patrickpfeifer.net",
    description="RNA-seq analysis | web frontend",
    license="Mixed",
    packages=find_packages(),
    namespace_packages = ['rnaseqlyze'],
    test_suite='nose.collector',
    include_package_data=True,
    zip_safe=False,
    data_files=[('', [
        'development.ini',
        'production.ini',
    ])],
    setup_requires=[
        'nose >= 1.1.2',
        'distribute >= 0.6.14',
        'setuptools_git >= 0.3',
    ],
    install_requires=[
        'pyramid',
        'pyramid_tm',
        'pyramid_debugtoolbar',
        'SQLAlchemy',
        'transaction',
        'PasteScript',
        'zope.sqlalchemy',
        'rna-seqlyze == ' + version,
    ],
    entry_points={
        'paste.app_factory': [
            'main = rnaseqlyze.web:main',
        ],
    },
)

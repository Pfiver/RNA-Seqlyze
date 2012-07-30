from setuptools import setup, find_packages

import rnaseqlyze
import rnaseqlyze.web

setup(
    name=rnaseqlyze.web.project_name,
    version=rnaseqlyze.__version__,
    author="Patrick Pfeifer",
    author_email="patrick@patrickpfeifer.net",
    description="RNA-seq analysis | web frontend",
    license="Mixed",
    packages=find_packages(),
    namespace_packages = ['rnaseqlyze'],
    test_suite='nose.collector',
    include_package_data=True,
    zip_safe=False,
    setup_requires=[
        'nose >= 1.1.2',
        'distribute >= 0.6.14',
        'setuptools_git >= 0.3',
        'pyramid',
        'pyramid_tm',
        'pyramid_debugtoolbar',
        'SQLAlchemy',
        'transaction',
        'PasteScript',
        'zope.sqlalchemy',
        'rna-seqlyze == ' + rnaseqlyze.__version__,
    ],
    entry_points={
        'paste.app_factory': [
            'main = rnaseqlyze.web:main',
        ],
    },
)

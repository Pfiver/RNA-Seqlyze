"""\
RNA-Seqlyze ApiDoc Generator

Usage:
    rnas-apidoc -h|--help
    rnas-apidoc [-s|--source] <path> ...
    
Generates one <package>.rst sphinx apidoc source file,
in the current directory, for each package found in <path>.

Options:
    -s --source    Use `literalinclude` in addition to `automodule`.
                   When using this option, the generated output will be
                   optimized for processing with the spinx latexpdf module that
                   generates a pdf document. Even without this option, when
                   using the html output module, the modules source code will
                   still be available in the generated Website, but not on the
                   same pages as the rest of the modules documentation
                   (controlled by the "html_show_sourcelink" option in
                   apidoc/conf.py).

"""
import os, sys
from pkgutil import walk_packages

def main():
    import docopt
    opts = docopt.docopt(__doc__)

    global pkg_tpl
    global mod_tpl
    if opts['--source']:
        pkg_tpl, mod_tpl = pkg_src_tpl, mod_src_tpl

    #: implicit args to write():
    #:  pkg.outfile, pgkpath, name, filename
    def write(tpl, **kwargs):
        pkg.outfile.write(tpl.format(
            name=name, path=pkgpath + os.sep + filename,
            equals="=" * len(name), dashes="-" * len(name), **kwargs))

    packages = {}
    Package = type('', (), {})

    # requirement/assumption:
    #  parent packages will come before their children
    for loader, name, is_pkg in walk_packages(opts['<path>']):

        pkgname = name.rsplit('.', 1)[0]
        pkgpath = os.path.relpath(loader.path, ".")

        if is_pkg:
            # note: 'pkgname' is actually the
            #       _parent_ package name in this case

            # associate non-root-packages with their parents
            if '.' in name:
                packages[pkgname].subpackages.append(name)

            # create & init a 'Package' object,
            # open <fully-qualified-package-name>.rst
            # and set filename to <package-name>/__init__.py
            pkg = Package()
            pkg.subpackages = []
            print "creating %s.rst" % name
            pkg.outfile = open(name + '.rst', 'w')
            filename = name.split('.')[-1] + os.sep + '__init__.py'

            # add the package to the list
            # and write the packages .rst file heading
            packages[name] = pkg
            write(pkg_tpl)

        else:
            # skip modules that are not part of any package, like setup.py
            if pkgname not in packages:
                continue

            # find the containing package and set filename to <module>.py
            pkg = packages[pkgname]
            filename = name[len(pkgname)+1:] + '.py'

            # append the doc entry for this module to the package .rst file
            write(mod_tpl)

    if not opts['--source']:
        for pkg in packages.values():
            if pkg.subpackages:
                write(sub_pkg_tpl, names="\n\t".join(pkg.subpackages))

pkg_tpl = """\
:mod:`{name}`
{equals}=======

.. automodule:: {name}

"""

pkg_src_tpl = """\
:mod:`{name}`
{equals}=======

:mod:`{name}`
{dashes}-------

.. automodule:: {name}

Source Code:

.. literalinclude:: {path}

"""

mod_tpl = """\
:mod:`{name}`
{dashes}-------

.. automodule:: {name}

"""

mod_src_tpl = """\
:mod:`{name}`
{dashes}-------

.. automodule:: {name}

Source Code:

.. literalinclude:: {path}

"""

sub_pkg_tpl = """\
Subpackages
-----------

.. toctree::
	:titlesonly:

	{names}

"""

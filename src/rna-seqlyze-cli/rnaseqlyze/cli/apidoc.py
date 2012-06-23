"""\
RNA-Seqlyze ApiDoc Generator

Usage:
    rnas-apidoc -h|--help
    rnas-apidoc [-s|--source] <path> ...
    
generates one <package>.rst sphinx apidoc source file,
in the current directory, for each package found in <path>

Options:
    -s --source    use `literalinclude` instead of `automodule` 
"""

def main():
    import os, sys
    from pkgutil import walk_packages

    import docopt
    opts = docopt.docopt(__doc__)

    Package = type('', (), {})
    packages = {}
    global mod_tpl
    if opts['--source']:
        mod_tpl = mod_src_tpl
    write = lambda tpl, **vals: \
            pkg.outfile.write(tpl.format(**vals))

    for loader, name, is_pkg in walk_packages(opts['<path>']):
        pkgname = name.rsplit('.', 1)[0]
        if is_pkg:
            if pkgname == name:
                modname, pkgname = pkgname, ''
            else:
                modname = name[len(pkgname)+1:]
                packages[pkgname].subpackages.append(name)
            pkg = packages[name] = Package()
            pkg.subpackages = []
            print "creating %s.rst" % name
            pkg.outfile = open(name + '.rst', 'w')
            write(pkg_tpl, name=name, equals="=" * len(name))
            filename = modname + os.sep + '__init__.py'

        else:
            # skip modules that are not part of any package, like setup.py
            if pkgname not in packages:
                continue
            pkg = packages[pkgname]
            filename = name[len(pkgname)+1:] + '.py'

        pkgpath = os.path.relpath(loader.path, ".")
        write(mod_tpl, name=name, dashes="-" * len(name),
                                  path=pkgpath + os.sep + filename)

    for pkg in packages.values():
        if pkg.subpackages:
            write(sub_pkg_tpl, names="\n\t".join(pkg.subpackages))

pkg_tpl = """\
Package `{name}`
{equals}==========

"""

mod_tpl = """\
:mod:`{name}`
{dashes}-------

.. automodule:: {name}
	:members:
	:undoc-members:
	:show-inheritance:

"""

mod_src_tpl = """\
:mod:`{name}`
{dashes}-------

.. literalinclude:: {path}

"""

sub_pkg_tpl = """\
Subpackages
-----------

.. toctree::

	{names}

"""

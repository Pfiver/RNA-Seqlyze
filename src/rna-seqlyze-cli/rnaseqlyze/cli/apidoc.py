import os, sys
from string import Template
from pkgutil import walk_packages

Package = type('Package', (object,), dict())

def main(argv=sys.argv):

    packages = {}
    module_tmpl = module
    subwrite = lambda tpl, **vals: \
            pkg.outfile.write(tpl.substitute(**vals))

    for loader, name, is_pkg in walk_packages(argv[1:]):
        pkgname = name.rsplit('.', 1)[0]
        if is_pkg:
            print name
            if pkgname == name:
                modname, pkgname = pkgname, ''
            else:
                modname = name[len(pkgname)+1:]
                packages[pkgname].subpackages.append(name)
            pkg = packages[name] = Package()
            pkg.subpackages = []
            pkg.outfile = open(name + '.rst', 'w')
            subwrite(package, name=name, equals="="*len(name))
            filename = modname + os.sep + '__init__.py'

        else:
            # skip modules that are not part of any package, like setup.py
            if pkgname not in packages:
                continue
            pkg = packages[pkgname]
            filename = name[len(pkgname)+1:] + '.py'

        if not pkgname:
            modpath, basepath = '', loader.path
        else:
            basepath = loader.path[:-len(pkgname)-1]
            modpath = loader.path[len(basepath)+1:] + os.sep

        pkgpath = os.path.basename(basepath) + os.sep
        subwrite(module_tmpl, name=name, dashes="-"*len(name),
                              path=pkgpath + modpath + filename)

    for pkg in packages.values():
        if pkg.subpackages:
            subwrite(subpackages, names="\n\t".join(pkg.subpackages))

package = Template("""\
Package `${name}`
${equals}==========

""")

module = Template("""\
:mod:`${name}`
${dashes}-------

.. automodule:: ${name}
	:members:
	:undoc-members:
	:show-inheritance:

""")

modulesrc = Template("""\
:mod:`${name}`
${dashes}-------

.. literalinclude:: ${path}

""")

subpackages = Template("""\
Subpackages
-----------

.. toctree::

	${names}

""")

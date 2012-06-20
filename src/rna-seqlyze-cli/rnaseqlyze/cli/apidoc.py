import sys
from string import Template
from pkgutil import walk_packages

Package = type('Package', (object,), dict())

def main(argv=sys.argv):
    packages = {}
    for loader, name, is_pkg in walk_packages(argv[1:]):
        pkgname = name.rsplit('.', 1)[0]
        if is_pkg:
            print name
            if pkgname != name:
                packages[pkgname].subpackages.append(name)
            pkg = packages[name] = Package()
            pkg.subpackages = []
            pkg.outfile = open(name + '.rst', 'w')
            pkg.outfile.write(package.substitute(name=name, equals="="*len(name)))
        else:
            # skip modules that are not part of any package
            if pkgname not in packages:
                continue
            pkg = packages[pkgname]
        pkg.outfile.write(module.substitute(name=name, dashes="-"*len(name)))
    for pkg in packages.values():
        if pkg.subpackages:
            pkg.outfile.write(subpackages.substitute(names="\n\t".join(pkg.subpackages)))

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

subpackages = Template("""\
Subpackages
-----------

.. toctree::

	${names}

""")

main()

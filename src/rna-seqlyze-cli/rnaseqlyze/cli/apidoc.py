import sys
from os import path
from string import Template
from pkgutil import walk_packages
from collections import defaultdict

Node = lambda: NodeType(Node)
NodeType = type('Node', (defaultdict,), dict())

def main(argv=sys.argv):
    tree = Node()
    tree.name = 'index'
    tree.is_pkg = True
    for loader, name, is_pkg in walk_packages(argv[1:]):
        if name == 'setup':
            continue
        node = tree
        for module in name.split('.'):
            node = node[module]
        node.is_pkg = is_pkg
        node.name = name
    walk(tree)

def walk(node):
    assert node.is_pkg == (len(node) > 0)
    if node.is_pkg:
        subs = []
        print("writing", node.name + '.rst')
        #outfile = open(node.name + '.rst', 'w')
        outfile = sys.stdout
        outfile.write(package.substitute(name=node.name, equals="="*len(node.name)))
        outfile.write(module.substitute(name=node.name, dashes="-"*len(node.name)))
        for k, v in node.items():
            if not v.is_pkg:
                outfile.write(module.substitute(name=v.name, dashes="-"*len(v.name)))
            else:
                subs.append(v.name)
                walk(v)
        if subs:
            outfile.write(subpackages.substitute(names="\n\t".join(subs)))
        #outfile.close()

package = Template("""\
Package `${name}`
${equals}========

""")

module = Template("""\
:mod:`${name}`
${dashes}-----

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

#!/usr/bin/python
# encoding: utf-8

data = {}

import re
import fileinput

for line in fileinput.input():

	line = line.strip()

	if len(line) > 0:

		words = re.split(r'\s+', line)
		item, deps = words[0], words[1:]

		data.setdefault(item, set()).update(set(deps))

fileinput.close()

def toposort(data):
	for k, v in data.items():
		v.discard(k)  # Ignore self dependencies
	extra_items_in_deps = reduce(set.union, data.values()) - set(data.keys())
	data.update({item: set() for item in extra_items_in_deps})
	while True:
		ordered = set(item for item, dep in data.items() if not dep)
		if not ordered:
			break
		yield ' '.join(sorted(ordered))
		data = {item: (dep - ordered) for item, dep in data.items()
				if item not in ordered}
	assert not data, "A cyclic dependency exists amongst %r" % data

print('\n'.join(toposort(data)))

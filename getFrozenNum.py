#!/usr/bin/python
#
# getFrozenNum coord
# outputs atom numbers that are frozen

import sys
import re

file = sys.argv[1]
frozenCoords = []
counter = 1

with open(file) as coordFile:
    for line in coordFile:
        if re.match('^\$.*', line):
            pass
        else:
            if len(line.split()) == 5:
                frozenCoords.append(counter)
            counter += 1


frozenCoords = ','.join(map(str, frozenCoords))
print frozenCoords

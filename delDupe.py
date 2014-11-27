#!/usr/bin/python

# Usage:
# delDupe.py inConstr
# deletes any duplicate lines in inConstr without otherwise reorganizing the file.
# e.g.:
# line1
# line2
# line1
# line3
# becomes:
# line1
# line2
# line3

import sys
from subprocess import call

inConstr = sys.argv[1]

f = open("tempInConstr", 'w')

model = 0

lines = []

with open(inConstr) as workingFile:
    for ii in workingFile:
        if ii in lines:
            pass
        else:
            lines.append(ii)
            f.write(ii)

f.close()

call('mv tempInConstr %s' % (inConstr), shell=True)



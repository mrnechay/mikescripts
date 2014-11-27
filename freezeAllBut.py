#!/usr/bin/python

# Usage:
# allFrozenBut 2,3,4,5 11
# spits 1,6,7,8,9,10,11

import sys

input = sys.argv[1]
numberOfAtoms = sys.argv[2]
freezeList = []

for ii in range(int(numberOfAtoms)):
    if str(ii + 1) not in input.split(","):
        freezeList.append(str(ii+1))
print ','.join(freezeList)



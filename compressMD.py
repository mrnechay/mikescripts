#!/usr/bin/python

# Usage:
# compessMD.py movie.pdb output.pdb 10
# compresses movie.pdb to output.pdb, only saving 1 every 10 frames

import sys

moviePDB = sys.argv[1]
outputF = sys.argv[2]
timeStep = sys.argv[3]

f = open(outputF, 'w')

model = 0

with open(moviePDB) as movieFile:
    for ii in movieFile:
        if model % int(timeStep) == 0:
            f.write(ii)
        if "ENDMDL" in ii:
            model += 1

f.close()
        

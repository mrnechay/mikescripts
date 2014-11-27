#!/usr/bin/python

# Usage:
# cat.py movie.pdb totalmovie.pdb output.pdb
# appends movie.pdb to totalmovie.pdb as output.pdb

import sys

moviePDB = sys.argv[1]
totalmovie = sys.argv[2]
outputF = sys.argv[3]

import sys

f = open(outputF, 'w')

with open(totalmovie) as movieFile:
    for ii in movieFile:
        f.write(ii)

with open(moviePDB) as movieFile:
    for ii in movieFile:
        f.write(ii)


f.close()
        

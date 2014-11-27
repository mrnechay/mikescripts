#!/usr/bin/python

# This converts the DMD placeholder Zn into the metal you actually want.

# As originally written by Tony, this script replaced metals based on their order in the new.pdb
# I ran into problems when a metal from outside the QM region (and thus not present in chop.xyz
# was ordered first in the new.pdb but not present in chop.xyz

# I wrote this to base substitution on position instead. DMD should not have moved any metals
# (otherwise, this code will unfortunatley break down...)

# Argument: xyz file to convert 
# (probably chop.xyz or subset.xyz)

# Additional notes: 
# The file new.pdb must be present.
# To generate a proper new.pdb file, the names of your desired atoms must be correct
# beginning from your initial.pdb file. Write the 1-2 letter element abbreviation in
# the atom name column (3rd column) of the pdb file. Example:
# HETATM   37 FE   WE1 B   4      98.764  96.268  93.628  1.00  0.00          Zn
#             ^^
# will replace Zn inthe chop.xyz with iron.

import re
import sys
import os

metals = []
atomserialnospaces = []
atomName = []
xCoordinate = []
yCoordinate = []
zCoordinate = []
elementSymbol = []

# load relevant metal atoms into memory (all metals should all be Zn at this point, so that is what the regex assumes)
with open('new.pdb') as pdbFile:
    for line in pdbFile:
        if re.match('^HETATM[ \t]+\d+[ \t]*([A-Za-z]+)\d*[ \t]+.*[ \t]+Zn\s*$', line):
            atomName.append(line[12:16].replace(" ", ""))
            xCoordinate.append(line[30:38].replace(" ", ""))
            yCoordinate.append(line[38:46].replace(" ", ""))
            zCoordinate.append(line[46:54].replace(" ", ""))
            elementSymbol.append(line[76:78].replace(" ", ""))

chop = sys.argv[1]
tmp = open('_tmp','w')
lineUpdated = False

with open(chop) as chopFile:
    for line in chopFile:
        splitLine = line.split()
        if len(splitLine) == 4: #it must be an atom line
            for idx, ii in enumerate(atomName):
                print splitLine[1]
                print xCoordinate[idx]
                print
                if splitLine[1] == xCoordinate[idx] and splitLine[2] == yCoordinate[idx] and splitLine[3] == zCoordinate[idx]:
                    print splitLine
                    lineUpdated = True
                    splitLine[0] = atomName[idx]
                    line2 = "%s      %s %s %s" % (splitLine[0], splitLine[1], splitLine[2], splitLine[3])
                    # ^ This will modify the line of the chopFile if it matched the positions of a known metal that needs to be replaced
        if lineUpdated:
            tmp.write(line2)
            lineUpdated = False
        else:
            tmp.write(line)

tmp.close()
os.rename('_tmp', chop)

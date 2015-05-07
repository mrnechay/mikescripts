#!/usr/bin/python

# didn't get very far. oh well.
import sys

file = sys.argv[1]

f = open('_temp.pdb')

with open(file + '.save') as saveFile:
    with open(file) as pdbFile:
        for saveLine in saveFile:
            for pdbLine in pdbFile:
                if "ATOM" in pdbLine.split(" ")[0] or "HETATM" in pdbLine.split(" ")[0]
                    if saveLine[21:26] == pdbLine[21:26] and saveLine[17:20] != pdbLine[17:20]:
                        if pdbLine
                    else:
                        pdbLine.write(line2)
                else:
                    

atom = []
atomserialnumber = []
atomserialnospaces = []
atomname = []
alternatelocationindicator = []
residuename = []
chainidentifier = []
residuesequencenumber = []
codeforinsertionofresidues = []
Xcoordinate = []
Ycoordinate = []
Zcoordinate = []
elementsymbol = []


COUNTER = 1
with open(pdbFile) as proteinfile:
    for ii in proteinfile:
        if "ATOM" in ii.split(" ")[0] or "HETATM" in ii.split(" ")[0]:
            atom.append(ii[0:6])
            atomserialnumber.append("%5i" % COUNTER)
            atomserialnospaces.append(int(COUNTER))
            atomname.append(ii[12:16])
            alternatelocationindicator.append(ii[16:17])
            residuename.append(ii[17:20])
            chainidentifier.append(ii[21:22])
            residuesequencenumber.append(ii[22:26])
            codeforinsertionofresidues.append(ii[26:27])
            Xcoordinate.append(float(ii[30:38]))
            Ycoordinate.append(float(ii[38:46]))
            Zcoordinate.append(float(ii[46:54]))
            elementsymbol.append(ii[76:78])

            atomsbonded = bondlist[COUNTER-1][1:]
            COUNTER += 1


        

#!/usr/bin/python

# This converts the DMD placeholder Zn into the metal you actually want.

# Argument: xyz file to convert 
# (probably chop.xyz or subset.xyz)

# Additional notes: 
# The file new.pdb must be present.
# To generate a proper new.pdb file, the names of your desired atoms must be correct
# beginning from your initial.pdb file. Write the 1-2 letter element abbreviation in
# the atom name column (3rd column) of the pdb file. Example:
# HETATM   37 Fe   WE1 B   4      98.764  96.268  93.628  1.00  0.00          Zn
#             ^^
# will replace Zn with iron.

import re
import sys
import os

metals = []
pattern = re.compile('^HETATM[ \t]+\d+[ \t]*([A-Za-z]+)\d*[ \t]+.*[ \t]+Zn\s*$')
f = open('new.pdb')
for line in f:
   match = re.match(pattern, line)
   if match is not None:
#      print match.group(1)
      metals.append(match.group(1))
f.close()
#print metals
i = 0
coordfile = sys.argv[1]
f = open(coordfile)
tmp = open('_tmp','w')
pattern2 = re.compile('[ \t]*[Zz][Nn]')
for line in f:
   if re.match(pattern2,line):
      line2 = re.sub(pattern2, metals[i], line)
      tmp.write(line2)
      i += 1
#      print line2
#      print i
   else:
      tmp.write(line)
f.close()
tmp.close()
os.rename('_tmp', coordfile)

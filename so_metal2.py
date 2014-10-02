#!/usr/bin/python

# This converts all metals back to Zn.

# Same caveats as with so_metal1.py apply.

# Rather than parsing new.pdb to determine what metals are present, this task could also be
# accomplished by matching all non-metals (anything not in [HCNOShcnos] and whatever else
# DMD recognizes) and replacing them with Zn. That sounds more elegant, but I am too lazy to
# look up the list of atoms that DMD recognizes.

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
for line in f:
   if re.match(str(metals[i]),line):
      line2 = re.sub(str(metals[i]), 'Zn', line)
      tmp.write(line2)
      i += 1
   else:
      tmp.write(line)
f.close()
tmp.close()
os.rename('_tmp', coordfile)

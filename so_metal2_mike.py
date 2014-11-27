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

chop = sys.argv[1]
tmp = open('_tmp','w')

with open(chop) as chopFile:
    for line in chopFile:
        splitLine = line.split()
        if len(splitLine) == 4: # we know it is an atom entry
            if splitLine[0].lower() in ["li", "be", "na", "mg", "k", "ca", "sc", "ti", "v", "cr", "mn", "fe", "co", "ni", "cu", "zn", "rb", "sr", "y", "zr", "nb", "mo", "tc", "ru", "pd", "ag", "cd", "cs", "ba", "hf", "ta", "w", "re", "os", "ir", "pt", "au", "hg", "la", "ce", "pr", "nd", "pm", "sm", "eu", "gd", "tb", "dy", "ho", "er", "tm", "yb", "lu"]:
                line2 = "Zn  %11s   %10s   %10s\n" % (splitLine[1], splitLine[2], splitLine[3])
                tmp.write(line2)
            else:
                tmp.write(line)
        else:
            tmp.write(line)
tmp.close()
os.rename('_tmp', chop)

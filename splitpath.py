#!/bin/bash

import re

"""
with open("path.xyz") as pathfile:
   for line in pathfile:
      numlines += 1
      if pat in line:
         break

numframes = sum(1 for line in open('path.xyz'))/numlines

print numframes

lastframe = 'f%d.pdb' % numframes
"""


def files():
   n = 0
   while True:
      n += 1
      yield open('number%d.xyz' % n, 'w') # This yield makes this a "generator function" which will continue with next() right where it left off

fs = files()
outfile = next(fs)
first = True

with open("path.xyz") as infile:
   for line in infile:
      isstartofxyz = re.search('^  +[0-9]', line)
      if not isstartofxyz:
         outfile.write(line)
      else:
         if first:
            first = False
            outfile.write(line)
         else:
            outfile.close()
            outfile = next(fs)
            outfile.write(line)

#!/usr/bin/python

import sys

filename = sys.argv[1]

HBXlist = []

with file(filename) as inputfile:
    for ii in inputfile:
        if ii[0:7] == '#frozen' and 'HBX' in ii:
            HBXlist.append(ii[8:17])

f = open('input2', 'w')

with file(filename) as inputfile:
    for ii in inputfile:
        if ii[0:7] == '#frozen' and 'HBX' not in ii and ii[8:17] in HBXlist:
            myList = list(ii.strip("\n"))
            myList.append(' HBX\n')
            jj = ''.join(myList)
            f.write(jj)
        elif ii[0:7] == '#frozen' and 'HBX' in ii:
            pass
        else:
            f.write(ii)

f.close()
        
        

#!/usr/bin/python

# usage:
# freezeatoms.py coord 1,2,5,10,20
# will freeze atoms 1, 2, 5, 10, and 20 in 'coord' file

import sys, os

coord = sys.argv[1]
frozenList = [int(x) for x in sys.argv[2].split(',')]

cartesianSection = False
atomIndex = 1

newCoord = open("__coord",'w')

with open(coord) as coordFile:
    for line in coordFile:
        if line == '$coord\n':
            cartesianSection = True
            newCoord.write(line)
        elif line[0] == '$':
            cartesianSection = False
            newCoord.write(line)
        elif cartesianSection == True:
            ls = filter(None, [x.strip().split(' ') for x in line.split('\n') if x.strip()][0])
            if atomIndex in frozenList:
                newCoord.write("%20.14f  %20.14f  %20.14f  %5s f\n" % (float(ls[0]), float(ls[1]), float(ls[2]), ls[3]))
            else:
                newCoord.write("%20.14f  %20.14f  %20.14f  %5s\n" % (float(ls[0]), float(ls[1]), float(ls[2]), ls[3]))
            atomIndex += 1
        elif cartesianSection == False:
            newCoord.write(line)

newCoord.close()
os.rename("__coord", "coord") 
    

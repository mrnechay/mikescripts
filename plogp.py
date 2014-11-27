#!/usr/bin/python

import numpy as np
from subprocess import call

numofhits = np.load("numofhits.npy")

#print len(numofhits)
numframes = 0
for x, value in np.ndenumerate(numofhits):
    
    numframes += value
    mmax = int(np.amax(numofhits))
    sumSoFar = 0
    mmax += 1

for ii in range(mmax):
    thing = len(numofhits[numofhits == ii]) / numframes
    if (thing > 0):
        sumSoFar += thing * np.log(thing)
                
print "p * log(p) = %f" % sumSoFar
call(' echo "p * log(p) = %f" >> plogp' % (sumSoFar), shell=True) 

#!/usr/bin/python

workingDir = sys.argv[1]
pdbFile = sys.argv[2]
qmdmdInput = "%s.qmdmd" % (workingDir)
startRead = False

with open(qmdmdInput) as qmdmdFile:
    for line in qmdmdFile:
        if 'END QM/DMD PARAMETERS' in line:
            startRead = True
        


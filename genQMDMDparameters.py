#!/bin/user/python

# This script reads from the .qmdmd input file and generates the files needed
# to interface TURBOMOLE and DMD including chopping the desired regions.

import numpy as np
from subprocess import call
import copy
import re
from operator import itemgetter

qmdmdInput = "default.qmdmd"
pdbFile = "try.pdb"

# initiate variables/dictionaries
clusterOptions = {
                 'slots': 4,
                 'time': 24,
                 'exclusive': True
}

charge = 0
multiplicity = 1

DMDoptions = {
             'Ti': 0.1,
             'Tf': 0.1,
             'Clustersize': 3,
             'Stepsannealing': 40,
             'Moviedt': 10,
             'Timesteps': 10000,
             'MaxIter': 40,
             'Equilibrate': 'true',
             'Discard': 'true',
             'Converge': 'true',
}

qmResidueList = []



# Sections of input file
cutSection = False

# following loop brings all options from the input file into memory
with open(qmdmdInput) as infile:
    for line in infile:
        if line[0] == "#" or not line or line.isspace():
            pass
        elif cutSection:
            qmResidueList.append(line.split())
        else:
            if "  slots" in line:
                clusterOptions['slots'] = line.split()[-1]
            elif "  time" in line:
                clusterOptions['time'] = line.split()[-1]
            elif "  exclusive" in line:
                clusterOptions['exclusive'] = line.split()[-1]
            elif "  charge" in line:
                charge = line.split()[-1]
            elif "  multiplicity" in line:
                multiplicity = line.split()[-1]
            elif "Ti" in line:
                DMDoptions['Ti'] = line.split()[-1]
            elif "Tf" in line:
                DMDoptions['Tf'] = line.split()[-1]
            elif "Cluster" in line:
                DMDoptions['Clustersize'] = line.split()[-1]
            elif "annealing" in line:
                DMDoptions['Stepsannealing'] = line.split()[-1]
            elif "dt" in line:
                DMDoptions['Moviedt'] = line.split()[-1]
            elif "DMD time steps" in line:
                DMDoptions['Timesteps'] = line.split()[-1]            
            elif "Max Iteration" in line:
                DMDoptions['MaxIter'] = line.split()[-1]            
            elif "Equilibrate" in line:
                DMDoptions['Equilibrate'] = line.split()[-1]       
            elif "Discard" in line:
                DMDoptions['Discard'] = line.split()[-1]       
            elif "Converge" in line:
                DMDoptions['Converge'] = line.split()[-1]       
            elif "QM Residues" in line:
                cutSection = True

babelCall = 'babel %s pdb.mol2 2> /dev/null' % (pdbFile)

call(babelCall, shell=True)
mol2file = "pdb.mol2"

def getBondList(file)
###############################################################################
# Let's create a bond table we can easily reference later to help in chopping
# While we could determine them based on pdb name, the format isn't consistent
# The mol2 file has a bond table, but it is in an inconvenient format. Fixing:

babelCall = 'babel %s pdb.mol2 2> /dev/null' % (pdbFile)

call(babelCall, shell=True)
mol2file = "pdb.mol2"

bondarea = False
bondtable = []


with open(mol2file) as bondfile:
    for ii in bondfile:
        if bondarea:
            # take columns 1 and 2 which are the two atoms bonded
            filtered = filter(bool, ii.strip('\n').split(" "))[1:3]
            changeInt = [ int(x) for x in filtered ]
            bondtable.append(changeInt)
            bondtable.append(changeInt[::-1])
        if "ATOM" == ii[9:13]:
            bondarea = False
        elif "BOND" == ii[9:13]:
            bondarea = True
    bondlist = sorted(bondtable, key=itemgetter(0))
    bondlistcopy = copy.deepcopy(bondlist)
    
counter = 0

#this loop merges entries with the same leading atom
for ii in bondlistcopy:
    lesst = 1
    if counter == 0:
        pass
    while bondlistcopy[counter][0] == bondlistcopy[counter-lesst][0]:
        lesst += 1
        if bondlistcopy[counter][0] != bondlistcopy[counter-lesst][0]:
            bondlist[counter-lesst+1].append(bondlistcopy[counter][1])
    counter += 1


for ii in range(len(bondlistcopy)):
    greatert = 1
    if ii + 2 < len(bondlistcopy):
        while bondlistcopy[ii][0] == bondlistcopy[ii+greatert][0]:
            bondlist[ii+greatert] = []
            if len(bondlistcopy) <= ii + greatert + 1:
                break
            else:
                greatert += 1

bondlist = filter(None, bondlist)

#add missing entries

goagain = True
while goagain:
    goagain = False
    counter = 0
    for ii in bondlist:
        if ii[0] == counter + 1:
            counter += 1
        elif ii[0] < counter + 1:
            del bondlist[counter]
            goagain = True
            break
        else:
            bondlist.insert(counter,[counter+1])
            goagain = True
            break

# Simply reference bondlist[atomserial-1][5] for mol2 atom identification. Also, bondlist[atomserial-1][2:5] for coords

# Bond table now in memory
##########################

# Let's now define elementlist which gives element name for each (atomserialnumber - 1)

elementlist = []
with open(pdbFile) as proteinfile:
    for row in proteinfile:
        if "ATOM" in row.split(" ")[0] or "HETATM" in row.split(" ")[0]:
            elementlist.append(row[76:78])


######################################################
# The following loops brings the pdb file into memory:

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

            #print bondlist[2015:]
            #print COUNTER
            atomsbonded = bondlist[COUNTER-1][1:]
            elementsbonded = [] #start the list of elements bonded to this one.

            for jj in atomsbonded:
                elementsbonded.append(elementlist[jj-1])
            COUNTER += 1

qmRes = [None] * 100
qmList = [None] * 100
qmChop = [None] * 100
qmFrez = [None] * 100

resCounter = 0
for ii in qmResidueList:
    if len(ii[0].split("-")) == 1:
    # If we land here, we are adding a single residue to the QM region (no backbone)
        aminoname = ii[0][0:3]
        aminonum = ii[0][3:]


        if aminoname.isupper():
        # UPPER CASE means the entire side chain should be included in qm region.
            qmList[resCounter] = []
            qmChop[resCounter] = ["CA", "CB"]
            qmFrez[resCounter] = ["CB"]
            for ii in atomserialnospaces:
                if re.match('%s.*%s ' % (aminoname, aminonum), # regex match the input with the pdb file
                            '%s A *%s ' % (residuename[ii-1], residuesequencenumber[ii-1]), re.IGNORECASE):
                    qmRes[resCounter] = '%s %s%s' % (residuename[ii-1], chainidentifier[ii-1], residuesequencenumber[ii-1])
                    # ^ I put this here so that it only adds a residue to the list if it actually exists in the protein
                    if atomname[ii-1].replace(" ", "") not in ["N", "O", "C", "HN", "CA"]:
                        qmList[resCounter].append(atomname[ii-1].replace(" ", ""))


        elif aminoname.islower():
        # _lower case_ means only the head of the side chain should be included in the qm region.
        # I will enforce that the side chain up until the first branch (e.g. ring) or hetero atom will be cut out
            qmList[resCounter] = [] #initialize list of atoms to include in QM region
            qmChop[resCounter] = []
            for ii in atomserialnospaces:
                if re.match('%s.*%s ' % (aminoname, aminonum), # here we just need to find the alpha carbon so we can bond-walk from there
                            '%s A *%s ' % (residuename[ii-1], residuesequencenumber[ii-1]), re.IGNORECASE) and atomname[ii-1] == " CA ":
                    qmRes[resCounter] = '%s %s%s' % (residuename[ii-1], chainidentifier[ii-1], residuesequencenumber[ii-1])
                    qmSkip = [" C  ", " O  ", " N  ", " HN "] # We already know we don't want these backbone atoms in our qm region
                    previousAtomName = " C  " # we will be trying to step along the bonds starting at CA, coming from C
                    atomDirections = bondlist[ii-1][1:] # list of atoms connected to this one
                    currentAtomName = " CA "
                    for jj in atomDirections:
                    # Here we are "stepping" to the next atom along the chain
                        if atomname[jj-1] not in qmSkip:
                        # If we "stepped" on an atom we already know we want to skip, we would need to go another direction
                        # If we made it here, it was a successful step.
                            currentAtom = jj
                            previousAtomName = currentAtomName
                            currentAtomName = atomname[jj-1]
                            keepWalking = True
                            while keepWalking:
                                if len(bondlist[currentAtom-1]) == 1:
                                    qmChop[resCounter] = [previousAtomName.replace(" ", ""), currentAtomName.replace(" ", "")] # This must be the end of the line. cut here.
                                    qmFrez[resCounter] = [currentAtomName.replace(" ", "")]
                                    qmSkip.append(previousAtomName)
                                    keepWalking = False
                                elif currentAtomName not in qmSkip and (len(bondlist[currentAtom-1][1:]) >= 3 or elementsymbol[currentAtom-1].replace(" ", "") != "C"):
                                    qmSkip.append(previousAtomName) #let's start the qmRegion here! We don't want to add this one to the atoms to skip
                                    qmChop[resCounter] = [previousAtomName.replace(" ", ""), currentAtomName.replace(" ", "")]
                                    qmFrez[resCounter] = [currentAtomName.replace(" ", "")]
                                    keepWalking = False
                                else: # We need to keep walking if we made it here!
                                    qmSkip.append(previousAtomName)
                                    for mm in bondlist[currentAtom-1][1:]:
                                        if atomname[mm-1] not in qmSkip:
                                            previousAtomName = currentAtomName
                                            currentAtom = mm
                                            currentAtomName = atomname[mm-1]

            for ii in atomserialnospaces:
                if re.match('%s.*%s' % (aminoname, aminonum), # regex match the input with the pdb file
                            '%s A *%s' % (residuename[ii-1], residuesequencenumber[ii-1]), re.IGNORECASE) and atomname[ii-1] not in qmSkip:
                    qmList[resCounter].append(atomname[ii-1].replace(" ", ""))
        resCounter += 1
     
    elif len(ii[0].split("-")) == 2:
    # This is a range of residues which means we will include the entire side chains and all or part of the backbone along the sequence
    # First, expand range to all residue names
        begRes = ii[0].split("-")[0]
        endRes = ii[0].split("-")[1]
        begResEnd = re.match('[A-Z]{3}[0-9]+', begRes, re.IGNORECASE).end()
        endResEnd = re.match('[A-Z]{3}[0-9]+', begRes, re.IGNORECASE).end()
        begAminoName = begRes[0:3]
        begAminoNum = begRes[3:begResEnd]
        endAminoName = endRes[0:3]
        endAminoNum = endRes[3:endResEnd]
        begResOption = begRes[begResEnd:]
        endResOption = endRes[endResEnd:]
 
        for ii in atomserialnospaces:
            qmList[resCounter] = []
            qmChop[resCounter] = []
            if re.match('%s.*%s ' % (begAminoName, begAminoNum), # regex match user input with the pdb file
                        '%s A *%s ' % (residuename[ii-1], residuesequencenumber[ii-1]), re.IGNORECASE) and '%s %s%s' % (residuename[ii-1], chainidentifier[ii-1], residuesequencenumber[ii-1]) not in qmRes:
                qmRes[resCounter] = '%s %s%s' % (residuename[ii-1], chainidentifier[ii-1], residuesequencenumber[ii-1])

                if begResOption == "":
                    qmSkip = [" C  ", " O  "]
                    qmChop[resCounter] = ["CA", "C"]
                    qmFrez[resCounter] = ["C"]
                    for kk in atomserialnospaces:
                        if re.match('%s.*%s' % (begAminoName, begAminoNum), # regex match the input with the pdb file
                                    '%s A *%s' % (residuename[kk-1], residuesequencenumber[kk-1]), re.IGNORECASE) and atomname[kk-1] not in qmSkip:
                            qmList[resCounter].append(atomname[kk-1].replace(" ", ""))

                elif begResOption == "N":
                    qmList[resCounter].append(" N  ")
                    qmChop[resCounter] = ["CA", "N"]
                    qmFrez[resCounter] = ["N"]

                keepGoing = True
                resCounter += 1
                currentline = ii + 1
                while keepGoing:
                    if re.match('%s.*%s ' % (endAminoName, endAminoNum), # regex match user input with the pdb file
                                '%s A *%s ' % (residuename[currentline-1], residuesequencenumber[currentline-1]), re.IGNORECASE) and '%s %s%s' % (residuename[currentline-1], chainidentifier[currentline-1], residuesequencenumber[currentline-1]) not in qmRes:
                    # This sees the end residue and stops the loop
                        qmRes[resCounter] = '%s %s%s' % (residuename[currentline-1], chainidentifier[currentline-1], residuesequencenumber[currentline-1])
                        qmList[resCounter] = []
                        keepGoing = False
                        
                        if endResOption == "N":
                            qmSkip = [" N  "]
                            qmChop[resCounter] = ["N", "CA"]
                            qmFrez[resCounter] = ["CA"]
 
                            for kk in atomserialnospaces:
                                if re.match('%s.*%s' % (endAminoName, BegAminoNum), # regex match the input with the pdb file
                                            '%s A *%s' % (residuename[kk-1], residuesequencenumber[kk-1]), re.IGNORECASE) and atomname[kk-1] not in qmSkip:
                                    qmList[resCounter].append(atomname[kk-1].replace(" ", ""))

                        elif begResOption == "":
                            qmList[resCounter].append([" C  ", " O  "])
                            qmChop[resCounter] = ["CA", "C"]
                            qmFrez[resCounter] = ["C"]
                        resCounter += 1

                    else:
                    # This part of the loop adds all atoms of everything within the range specified
                        currentResidue = '%s A%s' % (residuename[currentline-1], residuesequencenumber[currentline-1])
                        if currentResidue not in qmRes:
                            qmList[resCounter] = []
                            qmChop[resCounter] = ["Don't chop"]
                            qmFrez[resCounter] = ["Don't freeze"]
                            for kk in atomserialnospaces:
                                if re.match('%s.*%s ' % (residuename[currentline-1], residuesequencenumber[currentline-1]), # regex match the input with the pdb file
                                            '%s A *%s ' % (residuename[kk-1], residuesequencenumber[kk-1]), re.IGNORECASE):
                                    qmRes[resCounter] = currentResidue
                                    qmList[resCounter].append(atomname[kk-1].replace(" ", ""))
                            resCounter += 1
                        currentline += 1


# let's reorder that list of residues

qmList = filter(None, qmList)
qmRes = filter(None, qmRes)
qmChop = filter(None, qmChop)

third = []
for ii in qmRes:
    third.append(int(ii[5:].replace(" ", "")))
qmRes = [x for (y,x) in sorted(zip(third,qmRes), key=lambda pair: pair[0])]
qmList = [x for (y,x) in sorted(zip(third,qmList), key=lambda pair: pair[0])]
qmChop = [x for (y,x) in sorted(zip(third,qmChop), key=lambda pair: pair[0])]
qmFrez = [x for (y,x) in sorted(zip(third,qmFrez), key=lambda pair: pair[0])]

print "sorted lists"
print "qmRes"
print qmRes
print "qmList"
print qmList
print "qmChop"
print qmChop
print "qmFrez"
print qmFrez

print " "
print " "

# OK great, so we've added all of the atoms from the new.pdb into the qm region. but we're still missing nonpolar hydrogens.
# we will arrange a file in a format Chimera will be able to add H's for us




qmList = filter(None, qmList)
qmRes = filter(None, qmRes)
qmChop = filter(None, qmChop)

call('rm tempqm.pdb', shell=True)

extractQMcall = 'extract_qm.sh %s x.pdb >> qmdmdsetup.log' % (pdbFile)
call(extractQMcall, shell=True)
call('chimera_addh.sh x.pdb h.pdb >> qmdmdsetup.log 2>>qmdmdsetup.err', shell=True)




# this doesn't really work, for retarded reasons
#for idx, ii in enumerate(qmRes):
#    for jj in qmList[idx]:
#        grepString = 'egrep "%s +%s" %s >> tempqm.pdb' % (jj, qmRes[idx], pdbFile)
#        call(grepString, shell=True)
#    if qmChop[idx][0] == "Don't chop" or len(qmChop) > idx + 1 and qmChop[idx+1][0] == "Don't chop":
#        pass
#    else:
#        call('echo "TER" >> tempqm.pdb', shell=True)
#call('echo "END" >> tempqm.pdb', shell=True)

#call('chimera_addh.sh tempqm.pdb tempqmH.pdb', shell=True)

# At this point, we have tempqmH.pdb which has all the nonpolar H's we need for our qm region!
# now we need to add these to the lists we had before.


#for idx, ii in enumerate(qmRes):
#    with open('tempqmH.pdb') as Hfile:
#        for jj in Hfile:
#            if re.match('%s ' % (ii), jj[17:27]) and jj[12:16].replace(" ", "") not in qmList[idx]:
#                qmList[idx].append(jj[12:16].replace(" ", ""))




qmList = filter(None, qmList)
qmRes = filter(None, qmRes)
print qmRes
print qmList

#            splitLine = line.split()
#            if len(splitLine[0].split("-")) == 2:
#            # This line indicates two residues, these residues and their backbone
#            # will be included in the QM region

print qmResidueList

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
substrateCutList = []


# Sections of input file
cutSection = False
substrateCut = False

# following loop brings all options from the input file into memory
with open(qmdmdInput) as infile:
    for line in infile:
        if line[0] == "#" or not line or line.isspace():
            pass
        elif substrateCut:
            substrateCutList.append(line.strip().replace(" ", "").split('-'))
        elif cutSection:
            if "Custom " in line:
                substrateCut = True
            else:
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

def getBondList(fileToDo):
###############################################################################
# Let's create a bond table we can easily reference later to help in chopping
# While we could determine them based on pdb name, the format isn't consistent
# The mol2 file has a bond table, but it is in an inconvenient format. Fixing:

    babelCall = 'babel %s pdb.mol2 2> /dev/null' % (fileToDo)

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

    return bondlist
    # Simply reference bondlist[atomserial-1][5] for mol2 atom identification. Also, bondlist[atomserial-1][2:5] for coords

bondlist = getBondList(pdbFile)
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

            atomsbonded = bondlist[COUNTER-1][1:]
            COUNTER += 1

# pdb file now in memory
######################################################

# Initialize variables that carry the information for input lists:

qmRes = [None] * 100
qmList = [None] * 100
qmChop = [None] * 100
qmFrez = [None] * 100


# resCounter goes up +1 for every residue in the list. We will now cycle through those and add all of the atoms, their chopped and frozen components to the above lists:

resCounter = 0
for ii in qmResidueList:
    if len(ii[0].split("-")) == 1:
    # If we land here, we are adding a single residue to the QM region (no backbone)
        aminoname = ii[0][0:3]
        aminonum = ii[0][3:]
        option = ii[2].split(",")[0].replace(" ", "")

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

            if option == 'FrzAA':
                for jj in qmList[resCounter]:
                    if jj not in qmFrez[resCounter]:
                        qmFrez[resCounter].append(jj)

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
  
            if option == 'FrzAA':
                for jj in qmList[resCounter]:
                    if jj not in qmFrez[resCounter]:
                        qmFrez[resCounter].append(jj)

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
        option = ii[2].split(",")[0].replace(" ", "")
 
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
                if option  == 'FrzBB':
                    for jj in qmList[resCounter]:
                        if jj not in qmFrez[resCounter] and jj in ["N", "O", "C", "HN", "CA", "HA"]:
                            qmFrez[resCounter].append(jj)

                if option  == 'FrzAA':
                    for jj in qmList[resCounter]:
                        if jj not in qmFrez[resCounter]:
                            qmFrez[resCounter].append(jj)

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
                                if re.match('%s.*%s' % (endAminoName, endAminoNum), # regex match the input with the pdb file
                                            '%s A *%s' % (residuename[kk-1], residuesequencenumber[kk-1]), re.IGNORECASE) and atomname[kk-1] not in qmSkip:
                                    qmList[resCounter].append(atomname[kk-1].replace(" ", ""))

                        elif begResOption == "":
                            qmList[resCounter].append("C")
                            qmList[resCounter].append("O")
                            qmChop[resCounter] = ["CA", "C"]
                            qmFrez[resCounter] = ["C"]

                        if option  == 'FrzBB':
                            for jj in qmList[resCounter]:
                                if jj not in qmFrez[resCounter] and jj in ["N", "O", "C", "HN", "CA", "HA"]:
                                    qmFrez[resCounter].append(jj)

                        if option  == 'FrzAA':
                            for jj in qmList[resCounter]:
                                if jj not in qmFrez[resCounter]:
                                    qmFrez[resCounter].append(jj)

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

                            if option  == 'FrzBB':
                                for jj in qmList[resCounter]:
                                    qmFrez[resCounter] = []
                                    if jj not in qmFrez[resCounter] and jj in ["N", "O", "C", "HN", "CA", "HA"]:
                                        qmFrez[resCounter].append(jj)

                            if option  == 'FrzAA':
                                for jj in qmList[resCounter]:
                                    qmFrez[resCounter] = []
                                    if jj not in qmFrez[resCounter]:
                                        qmFrez[resCounter].append(jj)

                            resCounter += 1
                        currentline += 1

#############################################
# Let's add metals and substrate if there any
metalPresent = False # False until existence is proven
subPresent = False # False until existence is proven

for ii in atomserialnospaces:
    if atom[ii-1] == 'HETATM' and elementsymbol[ii-1].lower().replace(" ", "") in ["li", "be", "na", "mg", "k", "ca", "sc", "ti", "v", "cr", "mn", "fe", "co", "ni", "cu", "zn", "rb", "sr", "y", "zr", "nb", "mo", "tc", "ru", "pd", "ag", "cd", "cs", "ba", "hf", "ta", "w", "re", "os", "ir", "pt", "au", "hg", "la", "ce", "pr", "nd", "pm", "sm", "eu", "gd", "tb", "dy", "ho", "er", "tm", "yb", "lu"]:
        metalPresent = True
        if "%s %s%s" % (residuename[ii-1], chainidentifier[ii-1], residuesequencenumber[ii-1]) not in qmRes:
            qmList[resCounter] = [atomname[ii-1].replace(" ", "")]
            qmChop[resCounter] = ["Don't chop"]
            qmFrez[resCounter] = ["Don't freeze"]
            qmRes[resCounter] = "%s %s%s" % (residuename[ii-1], chainidentifier[ii-1], residuesequencenumber[ii-1])
            resCounter += 1
    elif atom[ii-1] == 'HETATM' and residuename[ii-1] == 'SUB':
        if "%s %s%s" % (residuename[ii-1], chainidentifier[ii-1], residuesequencenumber[ii-1]) not in qmRes:
            qmRes[resCounter] = "%s %s%s" % (residuename[ii-1], chainidentifier[ii-1], residuesequencenumber[ii-1]) 
            qmList[resCounter] = [] #initialize list of atoms to include in QM region
            # In the substrateCutList defined by the user, the first atom will be part of the qm region and the second will be cut out, or skipped from the qm region along with all atoms that come after.
            qmList[resCounter].append(list(x[0] for x in substrateCutList)[0])
            qmSkip = []
            qmSkip.append(list(x[1] for x in substrateCutList)[0])
            qmChop[resCounter] = substrateCutList
            qmFrez[resCounter] = list(x[0] for x in substrateCutList)
            qmSkip = []

            for idx, kk in enumerate(list(x[1] for x in substrateCutList)):
            # this loops through the atoms designated as the start of the DMD region
                for oo in atomserialnospaces:
                # We will loop through again to find atomserialnumber for this atom
                    if atomname[oo-1].replace(" ", "") == kk:
                        atomNUMBER = oo

                atomDirections = bondlist[atomNUMBER-1][1:] # list of atoms connected to this one
                keepWalking = True
                nextSteps = []
                             
                while keepWalking:
                    for ll in atomDirections:
                        if atomname[ll-1].replace(" ", "") not in qmList[resCounter] and atomname[ll-1].replace(" ", "") not in qmSkip:
                            qmSkip.append(atomname[ll-1])
                            for pp in bondlist[ll-1][1:]: #moving forward!
                                nextSteps.append(atomserialnospaces[pp-1])
                    if nextSteps == []:
                        keepWalking = False
                    else:
                        atomDirections = nextSteps
                        nextSteps = []
                for mm in atomserialnospaces:
                    if '%s %s' % (residuename[mm-1], chainidentifier[mm-1]) == 'SUB B' and atomname[mm-1].replace(" ", "") not in qmSkip:
                        qmList[resCounter].append(atomname[mm-1].replace(" ", ""))

            # If the list of cuts along the substrate is empty, let's include all atoms in the QM region and not freeze or chop anything
            if qmSkip == []:
                qmFrez[resCounter] = "Don't freeze"
                qmChop[resCounter] = "Don't chop"
                for pp in atomserialnospaces:
                    if "%s %s%s" % (residuename[pp-1], chainidentifier[pp-1], residuesequencenumber[pp-1]):
                        qmList[resCounter].append(atomname[pp-1].replace(" ", ""))

# Metal and/or substrates added
##############################################

# let's reorder that list of residues. If not in order, QM/DMD would naively attach the later-defined basis sets to the wrong atoms

qmList = filter(None, qmList)
qmRes = filter(None, qmRes)
qmChop = filter(None, qmChop)

third = []

for ii in qmRes:
    if ii.split()[1] == 'B':
        third.append(1000 - int(ii[5:].replace(" ", "")))
    else:
        third.append(int(ii[5:].replace(" ", "")))

qmRes = [x for (y,x) in sorted(zip(third,qmRes), key=lambda pair: pair[0])]
qmList = [x for (y,x) in sorted(zip(third,qmList), key=lambda pair: pair[0])]
qmChop = [x for (y,x) in sorted(zip(third,qmChop), key=lambda pair: pair[0])]
qmFrez = [x for (y,x) in sorted(zip(third,qmFrez), key=lambda pair: pair[0])]

# Sorting done
################

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

# OK great, so we've added all of the atoms from the new.pdb into the qm region. 
# But we're still missing nonpolar hydrogens.
# We will arrange a file in a format Chimera will be able to add H's for us

qmList = filter(None, qmList)
qmRes = filter(None, qmRes)
qmChop = filter(None, qmChop)

call('if [ -f tempqm.pdb ]; then rm tempqm.pdb; fi', shell=True)

# this section just to get reslist ready for Manuel's script to extract qm region
call('if [ -f input ]; then rm input; fi', shell=True)

myString = ""
for ii in qmRes:
    myString = myString + '"%s" ' % (ii)

myString = list(myString)
myString[-1] = ")"
myString.insert(0, "export reslist=(")
myString = "".join(myString)

f = open('input','w')
f.write(myString) # python will convert \n to os.linesep
f.close()

extractQMcall = 'extract_qm.sh %s x.pdb >> qmdmdsetup.log' % (pdbFile)
call(extractQMcall, shell=True)
call('chimera_addh.sh x.pdb h.pdb >> qmdmdsetup.log 2>>qmdmdsetup.err', shell=True)


# nonpolar hydrogens have been generated. Let's now add them to the qm region

bondListH = getBondList('h.pdb')

Hatom = []
Hatomserialnumber = []
Hatomserialnospaces = []
Hatomname = []
Halternatelocationindicator = []
Hresiduename = []
Hchainidentifier = []
Hresiduesequencenumber = []
Hcodeforinsertionofresidues = []
HXcoordinate = []
HYcoordinate = []
HZcoordinate = []
Helementsymbol = []

COUNTER = 1
with open('h.pdb') as Hfile:
    for ii in Hfile:
        if "ATOM" in ii.split(" ")[0] or "HETATM" in ii.split(" ")[0]:
            Hatom.append(ii[0:6])
            Hatomserialnumber.append("%5i" % COUNTER)
            Hatomserialnospaces.append(int(COUNTER))
            Hatomname.append(ii[12:16])
            Halternatelocationindicator.append(ii[16:17])
            Hresiduename.append(ii[17:20])
            Hchainidentifier.append(ii[21:22])
            Hresiduesequencenumber.append(ii[22:26])
            Hcodeforinsertionofresidues.append(ii[26:27])
            HXcoordinate.append(float(ii[30:38]))
            HYcoordinate.append(float(ii[38:46]))
            HZcoordinate.append(float(ii[46:54]))
            Helementsymbol.append(ii[76:78])

            COUNTER += 1
        


for idx, ii in enumerate(qmRes): #list of residues, also idx is the current index for qmList
    for idx2, jj in enumerate(Hresiduename): #idx2 is also the index of all of the H lists defined above
        for idx3, kk in enumerate(qmList[idx]): 
            if qmRes[idx] == "%s %s%s" % (Hresiduename[idx2], Hchainidentifier[idx2], Hresiduesequencenumber[idx2]) and kk == Hatomname[idx2].replace(" ", ""):
                # OK, at this point we've matched an atom of the qm region with a line in the h.pdb file
                # Now, we need to see if this atom is attached to an H in the h.pdb via the Hbondlist we made earlier.
                for ll in bondListH[idx2][1:]:
                    if Helementsymbol[ll-1].replace(" ", "") == "H" and Helementsymbol[idx2].replace(" ", "") not in ["N", "O", "S"]:
                    # Now append this H atom to the qm region of this residue
                        if Hatomname[ll-1].replace(" ", "") not in qmList[idx]:
                            qmList[idx].append(Hatomname[ll-1].replace(" ", ""))

# nonpolar hydrogens added to qm region
####################################################################


print qmRes
print "qmList:"
print qmList   
                       
# Great, all of the info for the input file should now be ready.

inputList = "#list "
inputChop = "#chop "
inputFrozen = "#frozen "

tempList = []
tempChop = []
tempFrozen = []

for idx, ii in enumerate(qmRes):

    inputList = "#list "
    inputChop = "#chop "
    inputFrozen = "#frozen "

    inputList = inputList + qmRes[idx]
    inputChop = inputChop + qmRes[idx] + ' '
    inputFrozen = inputFrozen + qmRes[idx] + ' '

    for jj in qmList[idx]:
        inputList = inputList + ' ' + jj
    tempList.append(inputList)

    if isinstance(qmChop[idx][0], str):
        if qmChop[idx][0] == "Don't chop":
            pass
        else:
            inputChop = inputChop + ' ' + qmChop[idx][1] + ' ' + qmChop[idx][0] + ' HBX 0.7052'
            tempChop.append(inputChop)
            inputFrozen = inputFrozen + ' ' + qmChop[idx][0] + ' HBX'
            tempFrozen.append(inputFrozen)

    elif isinstance(qmChop[idx][0], list):
        for jj in qmChop[idx]:
            inputChop = inputChop + ' ' + jj[1] + ' ' + jj[0] + ' HBX 0.7052'
            tempChop.append(inputChop)
            inputFrozen = inputFrozen + ' ' + jj[0] + ' HBX'
            tempFrozen.append(inputFrozen)
            inputChop = "#chop " 

call('echo "" >> input', shell=True)
call('echo " " >> input', shell=True)

for ii in tempList:
    call('echo "%s" >> input' % (ii), shell=True)
call('echo " " >> input', shell=True)

for ii in tempChop:
    call('echo "%s" >> input' % (ii), shell=True)
call('echo " " >> input', shell=True)

for ii in tempFrozen:
    call('echo "%s" >> input' % (ii), shell=True)
call('echo " " >> input', shell=True)

call('echo "export Ti=%s" >> input' % DMDoptions['Ti'], shell=True)
call('echo "export Tf=%s" >> input' % DMDoptions['Tf'], shell=True)
call('echo "export n_cluster=%s" >> input' % DMDoptions['Clustersize'], shell=True)
call('echo "export steps_annealing=%s" >> input' % DMDoptions['Stepsannealing'], shell=True)
call('echo "export f_movie_dt=%s" >> input' % DMDoptions['Moviedt'], shell=True)
call('echo "export MAX_TIME=%s" >> input' % DMDoptions['Timesteps'], shell=True)
call('echo "export Iterations=%s" >> input' % DMDoptions['MaxIter'], shell=True)
call('echo "export equilibrate=%s" >> input' % DMDoptions['Equilibrate'], shell=True)
call('echo "export discard=%s" >> input' % DMDoptions['Discard'], shell=True)
call('echo "export convergeDMD=%s" >> input' % DMDoptions['Converge'], shell=True)
call('echo "export dmd_cores=%s" >> input' % clusterOptions['slots'], shell=True)




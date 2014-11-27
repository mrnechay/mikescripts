#!/usr/bin/python

# This script reads from the .qmdmd input file and generates the files needed
# to interface TURBOMOLE and DMD including chopping the desired regions.

import numpy as np
from subprocess import call
import copy
import re
from operator import itemgetter
import sys
import json

qmdmdInput = "default.qmdmd"
pdbFile = "try.pdb"
workingDir = sys.argv[1]
pdbFile = sys.argv[2]
qmdmdInput = "%s.qmdmd" % (workingDir)

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
             'Dispmagnitude': '0.3'
}

allResFree = False

qmResidueList = []
substrateCutList = []
metalExcludeList = []
sidechainExcludeList = []
metalExclude = False
sidechainExclude = False

# Sections of input file
cutSection = False
substrateCut = False

# following loop brings all options from the input file into memory
with open(qmdmdInput) as infile:
    for line in infile:
        if 'END QM/DMD PARAMATERS' in line.upper():
            break
        if line[0] == "#" or not line or line.isspace():
            pass
        elif sidechainExclude:
            sidechainExcludeList.append(line.strip().replace(" ", ""))
        elif metalExclude:
            if "Exclude the following sidechain" in line:
                sidechainExclude = True
            else:
                metalExcludeList.append(line.strip().replace(" ", "").split('='))
        elif substrateCut:
            if "Exclude the following metals" in line:
                metalExclude = True
            else:
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
                charge = int(line.split()[-1])
            elif "  multiplicity" in line:
                multiplicity = int(line.split()[-1])
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
            elif "  Disp. magnitude" in line:
                DMDoptions['Dispmagnitude'] = float(line.split()[-1])
            elif "All residues free" in line:
                if line.split()[-1] == 'true':
                    allResFree = True
                else:
                    allResFree = False
            elif "QM Residues" in line:
                cutSection = True

# Let's first loop through the residue list to match protonation and deprotonation requests

# The previous script, setupqmdmd.mike, recognized incorrect HIS protonation and corrected it.
# If that was done, let's copy it over so the correction becomes the default protonation state.

allRes = []
allProtOpt = []
inConstrProtDeprot = []

#code refactored to not need below
#with open('inConstr') as inConstrFile:
#    for ii in inConstrFile:
#        if 'Protonate' in ii or 'Deprotonate' in ii:
#            print 'found protonation entry'
#            print ii
#            inConstrProtDeprot.append(ii)

# HIS correction has been loaded. On to user options

for idx, ii in enumerate(qmResidueList):
    if len(ii[-1].split("-")) > 1:
        allRes.append(ii[0].split("-")[0])
        allRes.append(ii[0].split("-")[0])
        allProtOpt.append(ii[-1].split("-")[0])
        allProtOpt.append(ii[-1].split("-")[1])
    else:
        allRes.append(ii[0])
        allProtOpt.append(ii[-1])

for idx, ii in enumerate(allRes):
    aminoname = ii[0][0:3].upper()
    aminonum = ii[0][3:]
    if allProtOpt[idx][0:4] == 'Prot':
        if aminoname == 'HIS':
            inConstrProtDeprot.append('Protonate A:%d:ND1' % (aminonum))
        if aminoname == 'ASP':
            inConstrProtDeprot.append('Protonate A:%d:OD%s' % (aminonum, allProtOpt[idx][4]))
        if aminoname == 'SER':
            inConstrProtDeprot.append('Protonate A:%d:OG' % (aminonum))
        if aminoname == 'THR':
            inConstrProtDeprot.append('Protonate A:%d:OG1' % (aminonum))
        if aminoname == 'ASP':
            if allProtOpt[idx][4] == 1:
                inConstrProtDeprot.append('Protonate A:%d:OD1' % (aminonum))
            elif allProtOpt[idx][4] == 2:
                inConstrProtDeprot.append('Protonate A:%d:OD2' % (aminonum))
        if aminoname == 'GLU':
            if allProtOpt[idx][4] == 1:
                inConstrProtDeprot.append('Protonate A:%d:OE1' % (aminonum))
            elif allProtOpt[idx][4] == 2:
                inConstrProtDeprot.append('Protonate A:%d:OE2' % (aminonum))
        if aminoname == 'ASN':
            if allProtOpt[idx][4] == 1:
                inConstrProtDeprot.append('Protonate A:%d:OD1' % (aminonum))
            elif allProtOpt[idx][4] == 2:
                inConstrProtDeprot.append('Protonate A:%d:ND2' % (aminonum))
        if aminoname == 'GLN':
            if allProtOpt[idx][4] == 1:
                inConstrProtDeprot.append('Protonate A:%d:OE1' % (aminonum))
            elif allProtOpt[idx][4] == 2:
                inConstrProtDeprot.append('Protonate A:%d:NE2' % (aminonum))
        if aminoname == 'TYR':
            inConstrProtDeprot.append('Protonate A:%d:OH' % (aminonum))
        if aminoname == 'TRP':
            inConstrProtDeprot.append('Protonate A:%d:NE1' % (aminonum))
    elif allProtOpt[idx][0:5] == 'DProt':
        if aminoname == 'ARG':
            inConstrProtDeprot.append('Deprotonate A:%d:NH1' % (aminonum))
        elif aminoname == 'HIS':
            inConstrProtDeprot.append('Deprotonate A:%d:NE1' % (aminonum))
        elif aminoname == 'LYS':
            inConstrProtDeprot.append('Deprotonate A:%d:NZ' % (aminonum))
        elif aminoname == 'CYS':
            inConstrProtDeprot.append('Deprotonate A:%d:SF' % (aminonum))
        elif aminoname == 'SEC':
            inConstrProtDeprot.append('Deprotonate A:%d:SEF' % (aminonum))
        elif aminoname == 'SER':
            inConstrProtDeprot.append('Deprotonate A:%d:OG' % (aminonum))
        elif aminoname == 'THR':
            inConstrProtDeprot.append('Deprotonate A:%d:OG1' % (aminonum))
        elif aminoname == 'ASN':
            inConstrProtDeprot.append('Deprotonate A:%d:ND2' % (aminonum))
        elif aminoname == 'GLN':
            inConstrProtDeprot.append('Deprotonate A:%d:NE2' % (aminonum))
        elif aminoname == 'TYR':
            inConstrProtDeprot.append('Deprotonate A:%d:OH' % (aminonum))


#I used to have the following lines to ensure a pdb file was DMD-compatible. BUT if following
# setupqmdmd.mike.new, it SHOULD be DMD-compatible at this point.
# Strangely, these lines give errors sometimes when an inConstr is present which is the step
# right before this one. i.e., even if pdb_to_pdbDMD.sh was successful at end of the first
# stage of setupqmdmd.mike, the generation of the inConstr file may break that script.
# in my test case, a zinc protein gave "atom overlap" errors for atom 5202 <-> 5202
# I have no idea why pdb_to_pdbDMD.sh was checking if the same atom overlaps with itself,
# but once the inConstr file is complete with Static terms for that metal, the problem goes
# away. very strange.
#call("pdb_to_pdbDMD.sh %s _temp.pdb >> qmdmdsetup.log 2>> qmdmdsetup.err" % (pdbFile), shell=True)
#call("cp _temp.pdb %s" % (pdbFile), shell=True)
#call("sed -i 's/Eh/ H/g' %s" % (pdbFile), shell=True)

####### Protonation and Deprotonation requests fulfilled ######################


def getBondList(fileToDo):
###############################################################################
# Let's create a bond table we can easily reference later to help in chopping
# While we could determine them based on pdb name, the format isn't consistent
# The mol2 file has a bond table, but it is in an inconvenient format. Fixing:

    babelCall = 'babel %s pdb.mol2 2> qmdmdsetup.err' % (fileToDo)
    
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

qmRes = [None] * 10000
qmList = [None] * 10000
qmChop = [None] * 10000
qmFrez = [None] * 10000


# resCounter goes up +1 for every residue in the list. We will now cycle through those and add all of the atoms, their chopped and frozen components to the above lists:

resCounter = 0
for ii in qmResidueList:
    if len(ii[0].split("-")) == 1:
    # If we land here, we are adding a single residue to the QM region (no backbone)
        aminoname = ii[0][0:3]
        aminonum = ii[0][3:]
        option = ii[2].split(",")[0].replace(" ", "")

        if aminoname.isupper() and '%s A%s' % (aminoname, aminonum) not in qmRes:
        # UPPER CASE means the entire side chain should be included in qm region.
            qmList[resCounter] = []
            qmChop[resCounter] = ["CA", "CB"]
            qmFrez[resCounter] = ["CB"]

            if aminoname == 'ARG':
                charge += 1
            if aminoname == 'LYS':
                charge += 1
            if aminoname == 'ASP':
                charge -= 1
            if aminoname == 'GLU':
                charge -= 1


            for jj in atomserialnospaces:
                if re.match('%s [A-Z] *%s ' % (aminoname, aminonum), # regex match the input with the pdb file
                            '%s %s%s ' % (residuename[jj-1], chainidentifier[jj-1], residuesequencenumber[jj-1]), re.IGNORECASE): 
                    qmRes[resCounter] = '%s %s%s' % (residuename[jj-1], chainidentifier[jj-1], residuesequencenumber[jj-1])
                    # ^ I put this here so that it only adds a residue to the list if it actually exists in the protein
                    if atomname[jj-1].replace(" ", "") not in ["N", "O", "C", "HN", "CA"]:
                        qmList[resCounter].append(atomname[jj-1].replace(" ", ""))

            if option == 'FrzAA':
                for jj in qmList[resCounter]:
                    if jj not in qmFrez[resCounter]:
                        qmFrez[resCounter].append(jj)

        elif aminoname.islower() and '%s A%s' % (aminoname, aminonum) not in qmRes:
        # _lower case_ means only the head of the side chain should be included in the qm region.
        # I will enforce that the side chain up until the first branch (e.g. ring) or hetero atom will be cut out
            qmList[resCounter] = [] #initialize list of atoms to include in QM region
            qmChop[resCounter] = []

            if aminoname == 'arg':
                charge += 1
            if aminoname == 'lys':
                charge += 1
            if aminoname == 'asp':
                charge -= 1
            if aminoname == 'glu':
                charge -= 1

            for ii in atomserialnospaces:
                if re.match('%s [A-Z] *%s ' % (aminoname, aminonum), # here we just need to find the alpha carbon so we can bond-walk from there
                            '%s %s%s ' % (residuename[ii-1], chainidentifier[ii-1], residuesequencenumber[ii-1]), re.IGNORECASE) and atomname[ii-1] == " CA ":

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
                if re.match('%s [A-Z] *%s' % (aminoname, aminonum), # regex match the input with the pdb file
                            '%s %s%s' % (residuename[ii-1], chainidentifier[ii-1], residuesequencenumber[ii-1]), re.IGNORECASE) and atomname[ii-1] not in qmSkip:
                    qmList[resCounter].append(atomname[ii-1].replace(" ", ""))
  
            if option == 'FrzAA':
                for jj in qmList[resCounter]:
                    if jj not in qmFrez[resCounter]:
                        qmFrez[resCounter].append(jj)

        resCounter += 1

    elif len(ii[0].split("-")) == 2:
        print ii
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
            if re.match('%s [A-Z] *%s ' % (begAminoName, begAminoNum), # regex match user input with the pdb file
                        '%s %s%s ' % (residuename[ii-1], chainidentifier[ii-1], residuesequencenumber[ii-1]), re.IGNORECASE) and '%s %s%s' % (residuename[ii-1], chainidentifier[ii-1], residuesequencenumber[ii-1]) not in qmRes:
                qmRes[resCounter] = '%s %s%s' % (residuename[ii-1], chainidentifier[ii-1], residuesequencenumber[ii-1])

                if begResOption == "N":
                    qmSkip = [" N  ", " HN "]
                    qmChop[resCounter] = ["N", "CA"]
                    qmFrez[resCounter] = ["CA"]
                    aminoname = begAminoName

                    if aminoname == 'ARG':
                        charge += 1
                    if aminoname == 'LYS':
                        charge += 1
                    if aminoname == 'ASP':
                        charge -= 1
                    if aminoname == 'GLU':
                        charge -= 1

                    for kk in atomserialnospaces:
                        if re.match('%s [A-Z] *%s' % (begAminoName, begAminoNum), # regex match the input with the pdb file
                                    '%s %s%s' % (residuename[kk-1], chainidentifier[kk-1], residuesequencenumber[kk-1]), re.IGNORECASE) and atomname[kk-1] not in qmSkip:
                            qmList[resCounter].append(atomname[kk-1].replace(" ", ""))

                elif begResOption == "":
                    qmList[resCounter].append("C")
                    qmList[resCounter].append("O")
                    qmChop[resCounter] = ["CA", "C"]
                    qmFrez[resCounter] = ["C"]

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

                call('echo "%s.*%s" >> tempthingy' % (endAminoName, endAminoNum), shell=True)
                while keepGoing:
                    if re.match('%s [A-Z] *%s ' % (endAminoName, endAminoNum), # regex match user input with the pdb file
                                '%s %s%s ' % (residuename[currentline-1], chainidentifier[currentline-1], residuesequencenumber[currentline-1]), re.IGNORECASE):
#Why did I have this appended to the above if statement?
#and '%s %s%s' % (residuename[currentline-1], chainidentifier[currentline-1], residuesequencenumber[currentline-1]) not in qmRes:
                    # This sees the end residue and stops the loop
                        qmRes[resCounter] = '%s %s%s' % (residuename[currentline-1], chainidentifier[currentline-1], residuesequencenumber[currentline-1])
                        qmList[resCounter] = []
                        keepGoing = False
                        call('echo "keepGoing is False" >> tempthingy', shell=True)
                        aminoname = endAminoName

                        if endResOption == "":
                            qmSkip = [" C  ", " O  "]
                            qmChop[resCounter] = ["C", "CA"]
                            qmFrez[resCounter] = ["CA"]

                            if aminoname == 'ARG':
                                charge += 1
                            if aminoname == 'LYS':
                                charge += 1
                            if aminoname == 'ASP':
                                charge -= 1
                            if aminoname == 'GLU':
                                charge -= 1

                            for kk in atomserialnospaces:
                                if re.match('%s [A-Z] *%s' % (endAminoName, endAminoNum), # regex match the input with the pdb file
                                            '%s %s%s' % (residuename[kk-1], chainidentifier[kk-1], residuesequencenumber[kk-1]), re.IGNORECASE) and atomname[kk-1] not in qmSkip:
                                    qmList[resCounter].append(atomname[kk-1].replace(" ", ""))

                        elif endResOption == "N":
                            qmList[resCounter].append("N")
                            qmList[resCounter].append("HN")
                            qmChop[resCounter] = ["CA", "N"]
                            qmFrez[resCounter] = ["N"]

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
                            if str("%s%s" % (residuename[currentline-1], residuesequencenumber[currentline-1])).replace(" ", "") in sidechainExcludeList:
                                qmList[resCounter] = ["N", "O", "C", "HN", "CA", "HA"]
                                qmChop[resCounter] = ["CB", "CA"]
                                qmFrez[resCounter] = ["CA"]
                                qmRes[resCounter] = currentResidue
                            else:
                                qmList[resCounter] = []
                                qmChop[resCounter] = ["Don't chop"]
                                qmFrez[resCounter] = ["Don't freeze"]
                                for kk in atomserialnospaces:
                                    if re.match('%s [A-Z] *%s ' % (residuename[currentline-1], residuesequencenumber[currentline-1]), # regex match the input with the pdb file
                                                '%s %s%s ' % (residuename[kk-1], chainidentifier[kk-1], residuesequencenumber[kk-1]), re.IGNORECASE):
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

                                aminoname = residuename[currentline-1]
                                if aminoname == 'ARG':
                                    charge += 1
                                if aminoname == 'LYS':
                                    charge += 1
                                if aminoname == 'ASP':
                                    charge -= 1
                                if aminoname == 'GLU':
                                    charge -= 1

                            resCounter += 1
                        currentline += 1

#############################################
# Let's add metals and substrate if there any
metalPresent = False # False until existence is proven
subPresent = False # False until existence is proven
metalList = []

for ii in atomserialnospaces:
    if atom[ii-1] == 'HETATM' and atomname[ii-1].lower().replace(" ", "") in ["li", "be", "na", "mg", "k", "ca", "sc", "ti", "v", "cr", "mn", "fe", "co", "ni", "cu", "zn", "rb", "sr", "y", "zr", "nb", "mo", "tc", "ru", "pd", "ag", "cd", "cs", "ba", "hf", "ta", "w", "re", "os", "ir", "pt", "au", "hg", "la", "ce", "pr", "nd", "pm", "sm", "eu", "gd", "tb", "dy", "ho", "er", "tm", "yb", "lu"]:
        if "%s %s%s" % (residuename[ii-1], chainidentifier[ii-1], residuesequencenumber[ii-1]) not in qmRes and residuename[ii-1].replace(" ", "") not in list(x[0] for x in metalExcludeList):
            metalPresent = True
            metalList.append(atomname[ii-1].replace(" ", ""))
            if atomname[ii-1].replace(" ", "") not in list(x[0] for x in metalExcludeList):
                qmList[resCounter] = [atomname[ii-1].replace(" ", "")]
                qmChop[resCounter] = ["Don't chop"]
                qmFrez[resCounter] = ["Don't freeze"]
                qmRes[resCounter] = "%s %s%s" % (residuename[ii-1], chainidentifier[ii-1], residuesequencenumber[ii-1])
                resCounter += 1
    elif atom[ii-1] == 'HETATM' and residuename[ii-1] == 'SUB':
        if "%s %s%s" % (residuename[ii-1], chainidentifier[ii-1], residuesequencenumber[ii-1]) not in qmRes and substrateCutList != []:
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
                    if '%s %s' % (residuename[mm-1], chainidentifier[mm-1]) == 'SUB B' and atomname[mm-1].replace(" ", "") not in qmSkip and atomname[mm-1].replace(" ", "") not in qmList[resCounter]:
                        qmList[resCounter].append(atomname[mm-1].replace(" ", ""))

        # If the list of cuts along the substrate is empty, let's include all atoms in the QM region and not freeze or chop anything
        elif substrateCutList == []:
            qmRes[resCounter] = "%s %s%s" % (residuename[ii-1], chainidentifier[ii-1], residuesequencenumber[ii-1])
            qmList[resCounter] = [] #initialize list of atoms to include in QM region
            qmFrez[resCounter] = ["Don't freeze"]
            qmChop[resCounter] = ["Don't chop"]
            for pp in atomserialnospaces:
                if "%s %s" % (residuename[pp-1], chainidentifier[pp-1]) == 'SUB B':
                    qmList[resCounter].append(atomname[pp-1].replace(" ", ""))

# Metal and/or substrates added
##############################################

# let's reorder that list of residues. If not in order, QM/DMD would naively attach the later-defined basis sets to the wrong atoms

qmList = filter(None, qmList)
qmRes = filter(None, qmRes)
qmChop = filter(None, qmChop)


third = [None] * len(qmRes)
newResOrder = []

currentCount = 1
for ii in atomserialnospaces:
    for idx, jj in enumerate(qmRes):
        if "%s %s%s" % (residuename[ii-1], chainidentifier[ii-1], residuesequencenumber[ii-1]) == jj and jj not in newResOrder:
            if "".join(re.split("[^a-zA-Z]*", residuename[ii-1])).lower() in ["li", "be", "na", "mg", "k", "ca", "sc", "ti", "v", "cr", "mn", "fe", "co", "ni", "cu", "zn", "rb", "sr", "y", "zr", "nb", "mo", "tc", "ru", "pd", "ag", "cd", "cs", "ba", "hf", "ta", "w", "re", "os", "ir", "pt", "au", "hg", "la", "ce", "pr", "nd", "pm", "sm", "eu", "gd", "tb", "dy", "ho", "er", "tm", "yb", "lu"]:
                third[idx] = currentCount + 200
                newResOrder.append(jj)
                currentCount += 1
            else:
                third[idx] = currentCount
                newResOrder.append(jj)
                currentCount += 1

qmRes = [x for (y,x) in sorted(zip(third,qmRes), key=lambda pair: pair[0])]
qmList = [x for (y,x) in sorted(zip(third,qmList), key=lambda pair: pair[0])]
qmChop = [x for (y,x) in sorted(zip(third,qmChop), key=lambda pair: pair[0])]
qmFrez = [x for (y,x) in sorted(zip(third,qmFrez), key=lambda pair: pair[0])]

# Sorting done
################

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
            inputChop = inputChop + ' ' + qmChop[idx][0] + ' ' + qmChop[idx][1] + ' HBX 0.7052'
            tempChop.append(inputChop)
            inputFrozen = inputFrozen + ' ' + qmChop[idx][1] + ' HBX'
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
call('echo "export dmd_cores=%s" >> input' % str(8), shell=True)
# Hardcoding "8" dmd cores above, because DMD does not parallelize well past 8 and the smallest nodes on hoffman have 8 cores

# input file generated
###########################################

###########################################
# Let's work on the inConstr file

inConstrAtompairrel = []
inConstrStatic = []

# Metal atoms need to be frozen
TERcount = 0
firstMetal = ""
with open(pdbFile) as pdbfileagain:
    for ii in pdbfileagain:
        if 'TER' in ii:
            TERcount += 1
        elif ii[12:14].lower().replace(" ", "") in ["li", "be", "na", "mg", "k", "ca", "sc", "ti", "v", "cr", "mn", "fe", "co", "ni", "cu", "zn", "rb", "sr", "y", "zr", "nb", "mo", "tc", "ru", "pd", "ag", "cd", "cs", "ba", "hf", "ta", "w", "re", "os", "ir", "pt", "au", "hg", "la", "ce", "pr", "nd", "pm", "sm", "eu", "gd", "tb", "dy", "ho", "er", "tm", "yb", "lu"]:
            inConstrStatic.append("Static %d.*.%s" % (TERcount + 1, ii[12:16].strip()))
            if firstMetal == "":
                firstMetal = "%d.%d.%s" % (TERcount + 1, int(ii[22:26].replace(" ", "")), ii[12:16].strip())

######### Metal atoms are now frozen in DMD ####################

# We need to know which atoms are in proximity with each of the metals. The script that is referenced will spit out those atoms closest. These be tagged as
# part of the "qm only" region, i.e. DMD is not allowed to move them (only QM)

# residues attached to metals
call("if [ -f metalsaminodist ]; then rm metalsaminodist; fi", shell=True)
for ii in metalList:
    call("atomstometal %s %s >> metalsaminodist" % (pdbFile, ii.upper()), shell=True)
resHeldStatic = []
with open('metalsaminodist') as metalDist:
    for ii in metalDist:
        atomStatic = ii.split(" ")[1]
        aminoname = ii.split(" ")[2][0:3]
        aminonum = ii.split(" ")[2][3:].strip()
        resHeldStatic.append("%s A%4d" % (aminoname, int(aminonum)))
        inConstrStatic.append("Static 1.%s.%s" % (aminonum, atomStatic))
        #I need to figure out which atom is connected to this
        #right now we have, e.g., atomStatic = OE1 and aminonum = 162
        for jj in atomserialnospaces:
            if atomname[jj-1].replace(" ", "") == atomStatic and int(residuesequencenumber[jj-1]) == int(aminonum):
                for kk in bondlist[jj-1][1:]:
                    if atomname[kk-1].lower().replace(" ", "") not in ["li", "be", "na", "mg", "k", "ca", "sc", "ti", "v", "cr", "mn", "fe", "co", "ni", "cu", "zn", "rb", "sr", "y", "zr", "nb", "mo", "tc", "ru", "pd", "ag", "cd", "cs", "ba", "hf", "ta", "w", "re", "os", "ir", "pt", "au", "hg", "la", "ce", "pr", "nd", "pm", "sm", "eu", "gd", "tb", "dy", "ho", "er", "tm", "yb", "lu"]:
                        inConstrAtompairrel.append("AtomPairRel 1.%s.%s %s -%f +%f" % (aminonum, atomname[kk-1].replace(" ", ""), firstMetal, str(float(DMDoptions['Dispmagnitude'])/5.0), str(float(DMDoptions['Dispmagnitude'])/5.0)))
        #Above is my attempt. let's see.

for idx, ii in enumerate(qmRes):
    if ii.split()[1] == 'A' and ii not in resHeldStatic:
    # We need to restrict the movement of this qm residue
        aminonum = ii[5:].replace(" ", "")
        for jj in qmList[idx]:
            if 'H' not in jj: #atom can't be nonpolar H, so to be easy let's just blanket this statement to all H
                inConstrAtompairrel.append("AtomPairRel 1.%s.%s %s -%f +%f" % (aminonum, jj, firstMetal, DMDoptions['Dispmagnitude'], DMDoptions['Dispmagnitude']))

# We should also, by default, leave the substrate as frozen to DMD. The only exception being if the user wanted to cut the substrate manually
TERcount = 0
with open(pdbFile) as pdbfileagain:
    for ii in pdbfileagain:
        if 'TER' in ii:
            TERcount += 1
        if 'SUB B' in ii:
            break
if substrateCutList == []:
    inConstrStatic.append("Static %d.*.*" % (TERcount + 1))
# if the user did cut the substrate, let's leave static all atoms in the list of qm atoms:
else:
    for idx, ii in enumerate(qmRes):
        if 'SUB B' in qmRes[idx]:
            for jj in qmList[idx]:
                inConstrStatic.append("Static %d.*.%s" % (TERcount + 1, jj)) 
########## Substrate atoms are now frozen in DMD ###############

# Let's print this stuff to the inConstr file!

call('echo " " >> inConstr', shell=True)

added1 = []
added2 = []
added3 = []
added4 = []

for ii in inConstrProtDeprot:
    if ii not in added1:
        added1.append(ii)
        call('echo "%s" >> inConstr' % (ii), shell=True)

call('echo " " >> inConstr', shell=True)
for ii in inConstrStatic:
    if ii not in added2:
        added2.append(ii)
        call('echo "%s" >> inConstr' % (ii), shell=True)

call('echo " " >> inConstr', shell=True)
for ii in inConstrAtompairrel:
    if ii not in added3:
        added3.append(ii)
        call('echo "%s" >> inConstr' % (ii), shell=True)

call('if [ -f qmRegionCharge ]; then rm qmRegionCharge; fi', shell=True)
call('echo "%d" >> qmRegionCharge' % (charge), shell=True)

qmResFile = open('temp/qmRes.txt', 'w')

qmListFile = open('temp/qmList.txt', 'w')

json.dump(qmRes, qmResFile)
json.dump(qmList, qmListFile)


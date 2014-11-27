#!/usr/bin/python

# Usage:

import json
from subprocess import call
from operator import itemgetter
import copy

qmResFile = open('temp/qmRes.txt', 'r')
qmListFile = open('temp/qmList.txt', 'r')

qmRes = json.load(qmResFile)
qmList = json.load(qmListFile)

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


tempList = []
tempChop = []
tempFrozen = []

for idx, ii in enumerate(qmRes):

    inputList = "#list "

    inputList = inputList + qmRes[idx]

    for jj in qmList[idx]:
        inputList = inputList + ' ' + jj
    tempList.append(inputList)

for ii in tempList:
    call('echo "%s" >> input2' % (ii), shell=True)


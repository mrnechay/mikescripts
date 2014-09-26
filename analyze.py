#!/usr/bin/python

import time
import os
import numpy as np
import scipy
import scipy.spatial
import copy
from subprocess import call

moviefile = "movie.pdb"

np.set_printoptions(threshold=np.nan)

Osp2residues = ["ASP", "GLU", "ASN", "GLN"]
Osp3residues = ["SER", "THR", "TYR"]
Narresidues = ["HIS"]

"""
definitions of groups used:

U - backbone hbond donors
V - backbone hbond acceptors
W - sidechain hbond acceptor, O atom sp2 hybridized (GLU, ASP)
X - sidechain hbond acceptor, O atom sp3 hybridized (SER)
Y - sidechain hbond acceptor, N atom aromatic (HIS)
ZO - sidechain hbond donor; H connected to O
ZN - sidechain hbond donor; H connected to N
"""

##################
# DMD Parameters!#

UV_OHdist_cutoff = [1.75, 2.50]
UV_NOdist_cutoff = [2.76, 3.50]
UV_HCdist_cutoff = [2.90, 3.75]

UW_OHdist_cutoff = [1.76, 3.00]
UW_CHdist_cutoff = [2.56, 4.25]
UW_ONdist_cutoff = [2.75, 4.00]

UX_OHdist_cutoff = [1.80, 3.00]
UX_CHdist_cutoff = [2.82, 4.50]
UX_ONdist_cutoff = [2.77, 4.00]

UY_NHdist_cutoff = [1.76, 3.00]
UY_CHdist_cutoff = [2.75, 4.25]
UY_NNdist_cutoff = [2.65, 4.00]

ZOV_OHdist_cutoff = [1.75, 3.00]
ZOV_OOdist_cutoff = [2.62, 4.00]
ZOV_HCdist_cutoff = [2.74, 4.25]

ZOW_OHdist_cutoff = [1.68, 3.00]
ZOW_OOdist_cutoff = [2.62, 4.00]
ZOW_HCdist_cutoff = [2.48, 4.25]

ZOX_OHdist_cutoff = [1.73, 3.00]
ZOX_OOdist_cutoff = [2.67, 4.00]
ZOX_HCdist_cutoff = [2.65, 4.50]

ZOY_NHdist_cutoff = [1.76, 3.00]
ZOY_CHdist_cutoff = [2.75, 4.25]
ZOY_NOdist_cutoff = [2.64, 4.00]

ZNV_OHdist_cutoff = [1.75, 3.00]
ZNV_NOdist_cutoff = [2.68, 4.00]
ZNV_HCdist_cutoff = [2.74, 4.25]

ZNW_OHdist_cutoff = [1.68, 3.00]
ZNW_ONdist_cutoff = [2.62, 4.00]
ZNW_HCdist_cutoff = [2.48, 4.25]

ZNX_OHdist_cutoff = [1.73, 3.00]
ZNX_ONdist_cutoff = [2.63, 4.00]
ZNX_HCdist_cutoff = [2.65, 4.50]

ZNY_NHdist_cutoff = [1.76, 3.00]
ZNY_CHdist_cutoff = [2.75, 4.25]
ZNY_NNdist_cutoff = [2.65, 4.00]

# End DMD parameters #
######################

#############################################
# Split movie.pdb into individual pdb files #

filescreated = []

numlines = 0
pat = 'ENDMDL'

with open("movie.pdb") as moviefile:
   for line in moviefile:
      numlines += 1
      if pat in line:
         break

numframes = sum(1 for line in open('movie.pdb'))/numlines

print numframes

lastframe = 'f%d.pdb' % numframes


def files():
   n = 0
   while True:
      n += 1
      filescreated.append('f%d.pdb' % n)
      yield open('f%d.pdb' % n, 'w') # This yield makes this a "generator function" which will continue with next() right where it left off

fs = files()
outfile = next(fs)
writeanother = False

first = True
with open("movie.pdb") as infile:
   for line in infile:
      if writeanother:
         outfile = next(fs)
         writeanother = False 
      if pat not in line:
         outfile.write(line)
      else:
         outfile.write(line)
         pdbnum = outfile.name
         pdbfile = pdbnum
         outfile.close()
         writeanother = True
         if first:
            call(["babel", outfile.name, "first.mol2"])
            first = False
            mol2file = "first.mol2"

# Movie being split to f1.pdb, f2.pdb, ...  #
# babel converted f1.pdb to first.mol2      #
#############################################

########################
# initialize variables #
            listofU = []
            listofV = []
            listofW = []
            listofX = []
            listofY = []
            listofZO = []
            listofZN = []

            listofUcoordsH = []
            listofVcoordsO = []
            listofWcoordsO = []
            listofXcoordsO = []
            listofYcoordsN = []
            listofZOcoordsH = []
            listofZNcoordsH = []

            listofUcoordsN = []
            listofVcoordsC = []
            listofWcoordsC = []
            listofXcoordsC = []
            listofYcoordsC = []
            listofZOcoordsO = []
            listofZNcoordsN = []


            elementlist = []
            bondtable = []
            bondarea = False


################################################################################################################
# Let's create a bond table we can _easily_ reference later to determine polar hydrogens.
# While we could determine them based on pdb name, the format isn't always consistent so there may be surprises.
# Plus, this would be the only way to determine them for the substrate.
# The mol2 file has a bond table, but it is in an inconvenient format. Let's fix that.

            with open(mol2file) as bondfile:
               for ii in bondfile:
                  if bondarea:
                     bondtable.append(filter(bool, ii.strip('\n').split(" "))[1:3])
                  if "ATOM" == ii[9:13]:
                     bondarea = False
                  elif "BOND" == ii[9:13]:
                     bondarea = True
               bondarray = np.array(bondtable)
               bondarray = bondarray.astype(np.int)
               secondbarray = bondarray[:,::-1]
               fullbarray = np.concatenate((bondarray,secondbarray), axis=0)

               sortedbonds = fullbarray[fullbarray[:,0].argsort()]

               bondlist = sortedbonds.tolist()
               bondlistcopy = copy.deepcopy(bondlist)

               counter = 0
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

# In case information about the residue is not available, we'll need to reference the mol2 atom identification to determine if O is sp2, sp3, or if N is ar

            triposatom = []
            checkforatom = []

            with open(mol2file) as bondfile:
               for jj in bondfile:
                  if checkforatom:
                     triposatom.append(filter(bool, jj.split(" ")))
                  if jj[9:13] == "ATOM":
                     checkforatom = True
                  elif jj[9:13] == "BOND":
                     checkforatom = False
                     del triposatom[-1]
            
# Simply reference triposatom[atomserial-1][5] for mol2 atom identification. Also, triposatom[atomserial-1][2:5] for coords

# Let's now define elementlist which gives element name for each (atomserialnumber - 1)

            with open(pdbfile) as proteinfile:
               for row in proteinfile:
                  if "ATOM" in row.split(" ")[0] or "HETATM" in row.split(" ")[0]:
                     elementlist.append(row[76:78])


# let's create our lists of atoms!
            COUNTER = 1
            with open(pdbfile) as proteinfile:
               for ii in proteinfile:
                  if "ATOM" in ii.split(" ")[0] or "HETATM" in ii.split(" ")[0]:
                     atom = ii[0:6]
                     atomserialnumber = "%5i" % COUNTER
                     atomserialnospaces = int(COUNTER)
                     atomname = ii[12:16]
                     alternatelocationindicator = ii[16:17]
                     residuename = ii[17:20]
                     chainidentifier = ii[21:22]
                     residuesequencenumber = ii[22:26]
                     codeforinsertionofresidues = ii[26:27]
                     Xcoordinate = float(ii[30:38])
                     Ycoordinate = float(ii[38:46])
                     Zcoordinate = float(ii[46:54])
                     elementsymbol = ii[76:78]
			
                     mol2identity = triposatom[atomserialnospaces-1][5]
			
                     if "ATOM" in atom and " HN " in atomname:
                        listofU.append(atomserialnospaces)

                     elif "ATOM" in atom and " O " in atomname:
                        listofV.append(atomserialnospaces)

                     elif (elementsymbol == " O" or elementsymbol == "O ") and (any(word in residuename for word in Osp2residues) or mol2identity == "O.2" or mol2identity == "O.co2"):
                        listofW.append(atomserialnospaces)

                     elif (elementsymbol == " O" or elementsymbol == "O ") and (any(word in residuename for word in Osp3residues) or mol2identity == "O.3"):
                        listofX.append(atomserialnospaces)

                     elif (elementsymbol == " N" or elementsymbol == "N ") and (any(word in residuename for word in Narresidues) or mol2identity == "N.ar"):
                        listofY.append(atomserialnospaces)

                     elif (elementsymbol == " H" or elementsymbol == "H "): 
                        atomsbonded = bondlist[atomserialnospaces-1][1:]
                        elementsbonded = [] #start the list of elements bonded to this one.
                        for jj in atomsbonded:
                           elementsbonded.append(elementlist[jj-1])
                        if (any(word in elementsbonded for word in ["O", " O", "O "])):
                           listofZO.append(atomserialnospaces)
                        elif (any(word in elementsbonded for word in ["N", " N", "N "])):
                           listofZN.append(atomserialnospaces)

                     COUNTER += 1

# This concludes a very important preliminary step - identifying the type of each hbond accepter/donor
# assuming no bonds break during the DMD simulaton, these types will stay the same throughout the run

# Now we need to iterate over all of the .pdb files we split up before.
# I'm going to repeat some code since there's a lot less we need to do now to loop over that again...

# intitialize all of the hbond arrays with zeros




            print "list of U"
            print listofU
            print "list of V"
            print listofV
            print "list of W"
            print listofW
            print "list of X"
            print listofX
            print "list of Y"
            print listofY
            print "list of ZO"
            print listofZO
            print "list of ZN"
            print listofZN




            UV_hbond_array = np.zeros((len(listofU),len(listofV)))
            UW_hbond_array = np.zeros((len(listofU),len(listofW)))
            UX_hbond_array = np.zeros((len(listofU),len(listofX)))
            UY_hbond_array = np.zeros((len(listofU),len(listofY)))
            ZOV_hbond_array = np.zeros((len(listofZO),len(listofV)))
            ZOW_hbond_array = np.zeros((len(listofZO),len(listofW)))
            ZOX_hbond_array = np.zeros((len(listofZO),len(listofX)))
            ZOY_hbond_array = np.zeros((len(listofZO),len(listofY)))
            ZNV_hbond_array = np.zeros((len(listofZN),len(listofV)))
            ZNW_hbond_array = np.zeros((len(listofZN),len(listofW)))
            ZNX_hbond_array = np.zeros((len(listofZN),len(listofX)))
            ZNY_hbond_array = np.zeros((len(listofZN),len(listofY)))
            
            numofhits = np.zeros(numframes)
            whole_3darray = []

         if pdbnum == lastframe:
            os.remove(lastframe)
            break #hack to avoid something weird

         UV_hbond_frame = np.zeros((len(listofU),len(listofV)))
         UW_hbond_frame = np.zeros((len(listofU),len(listofW)))
         UX_hbond_frame = np.zeros((len(listofU),len(listofX)))
         UY_hbond_frame = np.zeros((len(listofU),len(listofY)))
         ZOV_hbond_frame = np.zeros((len(listofZO),len(listofV)))
         ZOW_hbond_frame = np.zeros((len(listofZO),len(listofW)))
         ZOX_hbond_frame = np.zeros((len(listofZO),len(listofX)))
         ZOY_hbond_frame = np.zeros((len(listofZO),len(listofY)))
         ZNV_hbond_frame = np.zeros((len(listofZN),len(listofV)))
         ZNW_hbond_frame = np.zeros((len(listofZN),len(listofW)))
         ZNX_hbond_frame = np.zeros((len(listofZN),len(listofX)))
         ZNY_hbond_frame = np.zeros((len(listofZN),len(listofY)))

         listofUcoordsH = []
         listofVcoordsO = []
         listofWcoordsO = []
         listofXcoordsO = []
         listofYcoordsN = []
         listofZOcoordsH = []
         listofZNcoordsH = []

         listofUcoordsN = []
         listofVcoordsC = []
         listofWcoordsC = []
         listofXcoordsC = []
         listofXcoordsC_alt = []
         listofYcoordsC = []
         listofYcoordsC_alt = []
         listofZOcoordsO = []
         listofZNcoordsN = []

         pdbatom = []
         checkforatom = False

   #before we extracted coordinates from the mol2 since we needed the mol2 identifier as well, but now let's just use the pdb and avoid conversion to mol2
         with open(pdbnum) as pdbnumfile:
           for jj in pdbnumfile:
              if "ATOM" in jj.split(" ")[0] or "HETATM" in jj.split(" ")[0]:
                 pdbatom.append(filter(bool, jj[30:54].split(" ")))
   #pdb coordinates now cataloged in pdbatom[atomserial-1][0:3]


         for ii in listofU:
            for jj in bondlist[ii-1][1:]:
               if (any(word in elementlist[jj-1] for word in ["N ", " N"])):
                  listofUcoordsH.append(([float(kk) for kk in pdbatom[ii-1][0:3]]))
                  listofUcoordsN.append([float(kk) for kk in pdbatom[jj-1][0:3]])
                  break

         for ii in listofV:
            for jj in bondlist[ii-1][1:]:
               if (any(word in elementlist[jj-1] for word in ["C ", " C"])):
                  listofVcoordsO.append(([float(kk) for kk in pdbatom[ii-1][0:3]]))
                  listofVcoordsC.append([float(kk) for kk in pdbatom[jj-1][0:3]])
                  break

         for ii in listofW:
            for jj in bondlist[ii-1][1:]:
               if (any(word in elementlist[jj-1] for word in ["C ", " C"])):
                  listofWcoordsO.append(([float(kk) for kk in pdbatom[ii-1][0:3]]))
                  listofWcoordsC.append([float(kk) for kk in pdbatom[jj-1][0:3]])
                  break

         for ii in listofX:
            gonext = False
            wedone = False
      
            for jj in bondlist[ii-1][1:]:
               if wedone:
                  pass
               elif gonext:
                  if (any(word in elementlist[jj-1] for word in ["C ", " C"])):
                     listofXcoordsC_alt.append([float(kk) for kk in pdbatom[jj-1][0:3]])
                     wedone = True
                  elif (any(word in elementlist[jj-1] for word in ["H ", " H"])):
                     listofXcoordsC_alt.append([float(kk) for kk in pdbatom[jj-1][0:3]])
                     wedone = True
                  else:
                        listofXcoordsC_alt.append(listofXcoordsC[-1])
                        wedone = True
               elif (any(word in elementlist[jj-1] for word in ["C ", " C"])):
                  listofXcoordsO.append(([float(kk) for kk in pdbatom[ii-1][0:3]]))
                  listofXcoordsC.append([float(kk) for kk in pdbatom[jj-1][0:3]])
                  gonext = True
                  if jj == bondlist[ii-1][-1]:
                     listofXcoordsC_alt.append(listofXcoordsC[-1])
               elif (any(word in elementlist[jj-1] for word in ["H ", " H"])):
                  listofXcoordsO.append(([float(kk) for kk in pdbatom[ii-1][0:3]]))
                  listofXcoordsC.append([float(kk) for kk in pdbatom[jj-1][0:3]])
                  gonext = True
                  if jj == bondlist[ii-1][-1]:
                     listofXcoordsC_alt.append(listofXcoordsC[-1])
            if not wedone:
               listofXcoordsC_alt.append(listofXcoordsC[-1])
               wedone = True

         for ii in listofY:
            gonext = False
            wedone = False
            for jj in bondlist[ii-1][1:]:
               if wedone:
                  pass
               elif gonext:
                  if (any(word in elementlist[jj-1] for word in ["C ", " C"])):
                     listofYcoordsC_alt.append([float(kk) for kk in pdbatom[jj-1][0:3]])
                     wedone = True
                  else:
                     listofYcoordsC_alt.append(listofYcoordsC[-1])
                     wedone = True
               elif (any(word in elementlist[jj-1] for word in ["C ", " C"])):
                  listofYcoordsN.append(([float(kk) for kk in pdbatom[ii-1][0:3]]))
                  listofYcoordsC.append([float(kk) for kk in pdbatom[jj-1][0:3]])
                  gonext = True
                  if jj == bondlist[ii-1][-1]:
                     listofYcoordsC_alt.append([float(kk) for kk in pdbatom[jj-1][0:3]])
                     wedone = True
            if not wedone:
               listofYcoordsC_alt.append(listofYcoordsC[-1])

         for ii in listofZO:
            for jj in bondlist[ii-1][1:]:
               if (any(word in elementlist[jj-1] for word in ["O ", " O"])):
                  listofZOcoordsH.append(([float(kk) for kk in pdbatom[ii-1][0:3]]))
                  listofZOcoordsO.append([float(kk) for kk in pdbatom[jj-1][0:3]])
                  break


         for ii in listofZN:
            for jj in bondlist[ii-1][1:]:
               if (any(word in elementlist[jj-1] for word in ["N ", " N"])):
                  listofZNcoordsH.append(([float(kk) for kk in pdbatom[ii-1][0:3]]))
                  listofZNcoordsN.append([float(kk) for kk in pdbatom[jj-1][0:3]])
                  break

   # Calculate distances, courtesy of scipy's awesome "cdist" which gives a matrix of distances between every x with every y
   # Not all of the distances from these will be needed, e.g. an hbond check could "fail" after computing H---O or H---N distance,
   # But this cdist is the fastest distance calculator I could find and the result allows succinct code and logic from here

         if (len(listofUcoordsH) > 0 and len(listofVcoordsO) > 0):
            UV_OHdist = scipy.spatial.distance.cdist(listofUcoordsH,listofVcoordsO)
            UV_NOdist = scipy.spatial.distance.cdist(listofUcoordsN,listofVcoordsO)
            UV_HCdist = scipy.spatial.distance.cdist(listofUcoordsH,listofVcoordsC)

         if (len(listofUcoordsH) > 0 and len(listofWcoordsO) > 0):
            UW_OHdist = scipy.spatial.distance.cdist(listofUcoordsH,listofWcoordsO)
            UW_CHdist = scipy.spatial.distance.cdist(listofUcoordsH,listofWcoordsC)
            UW_ONdist = scipy.spatial.distance.cdist(listofUcoordsN,listofWcoordsO)

         if (len(listofUcoordsH) > 0 and len(listofXcoordsO) > 0):
            UX_OHdist = scipy.spatial.distance.cdist(listofUcoordsH,listofXcoordsO)
            UX_CHdist = scipy.spatial.distance.cdist(listofUcoordsH,listofXcoordsC)
            UX_CHdist_alt = scipy.spatial.distance.cdist(listofUcoordsH,listofXcoordsC_alt)
            UX_ONdist = scipy.spatial.distance.cdist(listofUcoordsN,listofXcoordsO)

         if (len(listofUcoordsH) > 0 and len(listofYcoordsN) > 0):
            UY_NHdist = scipy.spatial.distance.cdist(listofUcoordsH,listofYcoordsN)
            UY_CHdist = scipy.spatial.distance.cdist(listofUcoordsH,listofYcoordsC)
            UY_CHdist_alt = scipy.spatial.distance.cdist(listofUcoordsH,listofYcoordsC_alt)
            UY_NNdist = scipy.spatial.distance.cdist(listofUcoordsN,listofYcoordsN)
      
         if (len(listofZOcoordsH) > 0 and len(listofVcoordsO) > 0):
            ZOV_OHdist = scipy.spatial.distance.cdist(listofZOcoordsH,listofVcoordsO)
            ZOV_OOdist = scipy.spatial.distance.cdist(listofZOcoordsO,listofVcoordsO)
            ZOV_HCdist = scipy.spatial.distance.cdist(listofZOcoordsH,listofVcoordsC)

         if (len(listofZOcoordsH) > 0 and len(listofWcoordsO) > 0):
            ZOW_OHdist = scipy.spatial.distance.cdist(listofZOcoordsH,listofWcoordsO)
            ZOW_OOdist = scipy.spatial.distance.cdist(listofZOcoordsO,listofWcoordsO)
            ZOW_HCdist = scipy.spatial.distance.cdist(listofZOcoordsH,listofWcoordsC)

         if (len(listofZOcoordsH) > 0 and len(listofXcoordsO) > 0):
            ZOX_OHdist = scipy.spatial.distance.cdist(listofZOcoordsH,listofXcoordsO)
            ZOX_OOdist = scipy.spatial.distance.cdist(listofZOcoordsO,listofXcoordsO)
            ZOX_HCdist = scipy.spatial.distance.cdist(listofZOcoordsH,listofXcoordsC)
            ZOX_HCdist_alt = scipy.spatial.distance.cdist(listofZOcoordsH,listofXcoordsC_alt)

         if (len(listofZOcoordsH) > 0 and len(listofYcoordsN) > 0):
            ZOY_NHdist = scipy.spatial.distance.cdist(listofZOcoordsH,listofYcoordsN)
            ZOY_CHdist = scipy.spatial.distance.cdist(listofZOcoordsH,listofYcoordsC)
            ZOY_NOdist = scipy.spatial.distance.cdist(listofZOcoordsO,listofYcoordsN)
            ZOY_CHdist_alt = scipy.spatial.distance.cdist(listofZOcoordsH,listofYcoordsC_alt)

         if (len(listofZNcoordsH) > 0 and len(listofVcoordsO) > 0):
            ZNV_OHdist = scipy.spatial.distance.cdist(listofZNcoordsH,listofVcoordsO)
            ZNV_NOdist = scipy.spatial.distance.cdist(listofZNcoordsN,listofVcoordsO)
            ZNV_HCdist = scipy.spatial.distance.cdist(listofZNcoordsH,listofVcoordsC)

         if (len(listofZNcoordsH) > 0 and len(listofWcoordsO) > 0):
            ZNW_OHdist = scipy.spatial.distance.cdist(listofZNcoordsH,listofWcoordsO)
            ZNW_ONdist = scipy.spatial.distance.cdist(listofZNcoordsN,listofWcoordsO)
            ZNW_HCdist = scipy.spatial.distance.cdist(listofZNcoordsH,listofWcoordsC)

         if (len(listofZNcoordsH) > 0 and len(listofXcoordsO) > 0):
            ZNX_OHdist = scipy.spatial.distance.cdist(listofZNcoordsH,listofXcoordsO)
            ZNX_ONdist = scipy.spatial.distance.cdist(listofZNcoordsN,listofXcoordsO)
            ZNX_HCdist = scipy.spatial.distance.cdist(listofZNcoordsH,listofXcoordsC)
            ZNX_HCdist_alt = scipy.spatial.distance.cdist(listofZNcoordsH,listofXcoordsC_alt)

         if (len(listofZNcoordsH) > 0 and len(listofYcoordsN) > 0):
            ZNY_NHdist = scipy.spatial.distance.cdist(listofZNcoordsH,listofYcoordsN)
            ZNY_CHdist = scipy.spatial.distance.cdist(listofZNcoordsH,listofYcoordsC)
            ZNY_NNdist = scipy.spatial.distance.cdist(listofZNcoordsN,listofYcoordsN)
            ZNY_CHdist_alt = scipy.spatial.distance.cdist(listofZNcoordsH,listofYcoordsC_alt)

   
   # all pairs with atom pair distances within range get + 1

         if (len(listofUcoordsH) > 0 and len(listofVcoordsO) > 0):
            UV_hbond_array[(UV_OHdist > UV_OHdist_cutoff[0]) & (UV_OHdist < UV_OHdist_cutoff[1]) & \
                           (UV_NOdist > UV_NOdist_cutoff[0]) & (UV_NOdist < UV_NOdist_cutoff[1]) & \
                           (UV_HCdist > UV_HCdist_cutoff[0]) & (UV_HCdist < UV_HCdist_cutoff[1])] += 1

         if (len(listofUcoordsH) > 0 and len(listofWcoordsO) > 0):
            UW_hbond_array[(UW_OHdist > UW_OHdist_cutoff[0]) & (UW_OHdist < UW_OHdist_cutoff[1]) & \
                           (UW_CHdist > UW_CHdist_cutoff[0]) & (UW_CHdist < UW_CHdist_cutoff[1]) & \
                           (UW_ONdist > UW_ONdist_cutoff[0]) & (UW_ONdist < UW_ONdist_cutoff[1])] += 1
   
         if (len(listofUcoordsH) > 0 and len(listofXcoordsO) > 0):
            UX_hbond_array[(UX_OHdist > UX_OHdist_cutoff[0]) & (UX_OHdist < UX_OHdist_cutoff[1]) & \
                         (((UX_CHdist > UX_CHdist_cutoff[0]) & (UX_CHdist < UX_CHdist_cutoff[1])) | ((UX_CHdist_alt > UX_CHdist_cutoff[0]) & (UX_CHdist_alt < UX_CHdist_cutoff[1]))) & \
                           (UX_ONdist > UX_ONdist_cutoff[0]) & (UX_ONdist < UX_ONdist_cutoff[1])] += 1

         if (len(listofUcoordsH) > 0 and len(listofYcoordsN) > 0):
            UY_hbond_array[(UY_NHdist > UY_NHdist_cutoff[0]) & (UY_NHdist < UY_NHdist_cutoff[1]) & \
                         (((UY_CHdist > UY_CHdist_cutoff[0]) & (UY_CHdist < UY_CHdist_cutoff[1])) | ((UY_CHdist_alt > UY_CHdist_cutoff[0]) & (UY_CHdist_alt < UY_CHdist_cutoff[1]))) & \
                           (UY_NNdist > UY_NNdist_cutoff[0]) & (UY_NNdist < UY_NNdist_cutoff[1])] += 1
   
         if (len(listofZOcoordsH) > 0 and len(listofVcoordsO) > 0):
            ZOV_hbond_array[(ZOV_OHdist > ZOV_OHdist_cutoff[0]) & (ZOV_OHdist < ZOV_OHdist_cutoff[1]) & \
                            (ZOV_OOdist > ZOV_OOdist_cutoff[0]) & (ZOV_OOdist < ZOV_OOdist_cutoff[1]) & \
                            (ZOV_HCdist > ZOV_HCdist_cutoff[0]) & (ZOV_HCdist < ZOV_HCdist_cutoff[1])] += 1
   
         if (len(listofZOcoordsH) > 0 and len(listofWcoordsO) > 0):
            ZOW_hbond_array[(ZOW_OHdist > ZOW_OHdist_cutoff[0]) & (ZOW_OHdist < ZOW_OHdist_cutoff[1]) & \
                            (ZOW_OOdist > ZOW_OOdist_cutoff[0]) & (ZOW_OOdist < ZOW_OOdist_cutoff[1]) & \
                            (ZOW_HCdist > ZOW_HCdist_cutoff[0]) & (ZOW_HCdist < ZOW_HCdist_cutoff[1])] += 1
   
         if (len(listofZOcoordsH) > 0 and len(listofXcoordsO) > 0):   
            ZOX_hbond_array[(ZOX_OHdist > ZOX_OHdist_cutoff[0]) & (ZOX_OHdist < ZOX_OHdist_cutoff[1]) & \
                            (ZOX_OOdist > ZOX_OOdist_cutoff[0]) & (ZOX_OOdist < ZOX_OOdist_cutoff[1]) & \
                          (((ZOX_HCdist > ZOX_HCdist_cutoff[0]) & (ZOX_HCdist < ZOX_HCdist_cutoff[1])) | ((ZOX_HCdist_alt > ZOX_HCdist_cutoff[0]) & (ZOX_HCdist_alt < ZOX_HCdist_cutoff[1])))] += 1
   
         if (len(listofZOcoordsH) > 0 and len(listofYcoordsN) > 0):
            ZOY_hbond_array[(ZOY_NHdist > ZOY_NHdist_cutoff[0]) & (ZOY_NHdist < ZOY_NHdist_cutoff[1]) & \
                          (((ZOY_CHdist > ZOY_CHdist_cutoff[0]) & (ZOY_CHdist < ZOY_CHdist_cutoff[1])) | ((ZOY_CHdist_alt > ZOY_CHdist_cutoff[0]) & (ZOY_CHdist_alt < ZOY_CHdist_cutoff[1]))) & \
                            (ZOY_NOdist > ZOY_NOdist_cutoff[0]) & (ZOY_NOdist < ZOY_NOdist_cutoff[1])] += 1
   
         if (len(listofZNcoordsH) > 0 and len(listofVcoordsO) > 0):
            ZNV_hbond_array[(ZNV_OHdist > ZNV_OHdist_cutoff[0]) & (ZNV_OHdist < ZNV_OHdist_cutoff[1]) & \
                            (ZNV_NOdist > ZNV_NOdist_cutoff[0]) & (ZNV_NOdist < ZNV_NOdist_cutoff[1]) & \
                            (ZNV_HCdist > ZNV_HCdist_cutoff[0]) & (ZNV_HCdist < ZNV_HCdist_cutoff[1])] += 1
   
         if (len(listofZNcoordsH) > 0 and len(listofWcoordsO) > 0):
            ZNW_hbond_array[(ZNW_OHdist > ZNW_OHdist_cutoff[0]) & (ZNW_OHdist < ZNW_OHdist_cutoff[1]) & \
                            (ZNW_ONdist > ZNW_ONdist_cutoff[0]) & (ZNW_ONdist < ZNW_ONdist_cutoff[1]) & \
                            (ZNW_HCdist > ZNW_HCdist_cutoff[0]) & (ZNW_HCdist < ZNW_HCdist_cutoff[1])] += 1
   
         if (len(listofZNcoordsH) > 0 and len(listofXcoordsO) > 0):
            ZNX_hbond_array[(ZNX_OHdist > ZNX_OHdist_cutoff[0]) & (ZNX_OHdist < ZNX_OHdist_cutoff[1]) & \
                            (ZNX_ONdist > ZNX_ONdist_cutoff[0]) & (ZNX_ONdist < ZNX_ONdist_cutoff[1]) & \
                          (((ZNX_HCdist > ZNX_HCdist_cutoff[0]) & (ZNX_HCdist < ZNX_HCdist_cutoff[1])) | ((ZNX_HCdist > ZNX_HCdist_cutoff[0]) & (ZNX_HCdist < ZNX_HCdist_cutoff[1])))] += 1
   
         if (len(listofZNcoordsH) > 0 and len(listofYcoordsN) > 0):
            ZNY_hbond_array[(ZNY_NHdist > ZNY_NHdist_cutoff[0]) & (ZNY_NHdist < ZNY_NHdist_cutoff[1]) & \
                            (ZNY_CHdist > ZNY_CHdist_cutoff[0]) & (ZNY_CHdist < ZNY_CHdist_cutoff[1]) & \
                          (((ZNY_NNdist > ZNY_NNdist_cutoff[0]) & (ZNY_NNdist < ZNY_NNdist_cutoff[1])) | ((ZNY_NNdist > ZNY_NNdist_cutoff[0]) & (ZNY_NNdist < ZNY_NNdist_cutoff[1])))] += 1

   # 3d array, each slice is a frame
   
         if (len(listofUcoordsH) > 0 and len(listofVcoordsO) > 0):
            UV_hbond_frame[(UV_OHdist > UV_OHdist_cutoff[0]) & (UV_OHdist < UV_OHdist_cutoff[1]) & \
                           (UV_NOdist > UV_NOdist_cutoff[0]) & (UV_NOdist < UV_NOdist_cutoff[1]) & \
                           (UV_HCdist > UV_HCdist_cutoff[0]) & (UV_HCdist < UV_HCdist_cutoff[1])] = 1
   
         if (len(listofUcoordsH) > 0 and len(listofWcoordsO) > 0):
            UW_hbond_frame[(UW_OHdist > UW_OHdist_cutoff[0]) & (UW_OHdist < UW_OHdist_cutoff[1]) & \
                           (UW_CHdist > UW_CHdist_cutoff[0]) & (UW_CHdist < UW_CHdist_cutoff[1]) & \
                           (UW_ONdist > UW_ONdist_cutoff[0]) & (UW_ONdist < UW_ONdist_cutoff[1])] = 1
   
         if (len(listofUcoordsH) > 0 and len(listofXcoordsO) > 0):
            UX_hbond_frame[(UX_OHdist > UX_OHdist_cutoff[0]) & (UX_OHdist < UX_OHdist_cutoff[1]) & \
                         (((UX_CHdist > UX_CHdist_cutoff[0]) & (UX_CHdist < UX_CHdist_cutoff[1])) | ((UX_CHdist_alt > UX_CHdist_cutoff[0]) & (UX_CHdist_alt < UX_CHdist_cutoff[1]))) & \
                           (UX_ONdist > UX_ONdist_cutoff[0]) & (UX_ONdist < UX_ONdist_cutoff[1])] = 1
   
         if (len(listofUcoordsH) > 0 and len(listofYcoordsN) > 0):
            UY_hbond_frame[(UY_NHdist > UY_NHdist_cutoff[0]) & (UY_NHdist < UY_NHdist_cutoff[1]) & \
                         (((UY_CHdist > UY_CHdist_cutoff[0]) & (UY_CHdist < UY_CHdist_cutoff[1])) | ((UY_CHdist_alt > UY_CHdist_cutoff[0]) & (UY_CHdist_alt < UY_CHdist_cutoff[1]))) & \
                           (UY_NNdist > UY_NNdist_cutoff[0]) & (UY_NNdist < UY_NNdist_cutoff[1])] = 1
   
         if (len(listofZOcoordsH) > 0 and len(listofVcoordsO) > 0):
            ZOV_hbond_frame[(ZOV_OHdist > ZOV_OHdist_cutoff[0]) & (ZOV_OHdist < ZOV_OHdist_cutoff[1]) & \
                            (ZOV_OOdist > ZOV_OOdist_cutoff[0]) & (ZOV_OOdist < ZOV_OOdist_cutoff[1]) & \
                            (ZOV_HCdist > ZOV_HCdist_cutoff[0]) & (ZOV_HCdist < ZOV_HCdist_cutoff[1])] = 1
   
         if (len(listofZOcoordsH) > 0 and len(listofWcoordsO) > 0):
            ZOW_hbond_frame[(ZOW_OHdist > ZOW_OHdist_cutoff[0]) & (ZOW_OHdist < ZOW_OHdist_cutoff[1]) & \
                            (ZOW_OOdist > ZOW_OOdist_cutoff[0]) & (ZOW_OOdist < ZOW_OOdist_cutoff[1]) & \
                            (ZOW_HCdist > ZOW_HCdist_cutoff[0]) & (ZOW_HCdist < ZOW_HCdist_cutoff[1])] = 1
   
         if (len(listofZOcoordsH) > 0 and len(listofXcoordsO) > 0):
            ZOX_hbond_frame[(ZOX_OHdist > ZOX_OHdist_cutoff[0]) & (ZOX_OHdist < ZOX_OHdist_cutoff[1]) & \
                            (ZOX_OOdist > ZOX_OOdist_cutoff[0]) & (ZOX_OOdist < ZOX_OOdist_cutoff[1]) & \
                          (((ZOX_HCdist > ZOX_HCdist_cutoff[0]) & (ZOX_HCdist < ZOX_HCdist_cutoff[1])) | ((ZOX_HCdist_alt > ZOX_HCdist_cutoff[0]) & (ZOX_HCdist_alt < ZOX_HCdist_cutoff[1])))] = 1
   
         if (len(listofZOcoordsH) > 0 and len(listofYcoordsN) > 0):
            ZOY_hbond_frame[(ZOY_NHdist > ZOY_NHdist_cutoff[0]) & (ZOY_NHdist < ZOY_NHdist_cutoff[1]) & \
                          (((ZOY_CHdist > ZOY_CHdist_cutoff[0]) & (ZOY_CHdist < ZOY_CHdist_cutoff[1])) | ((ZOY_CHdist_alt > ZOY_CHdist_cutoff[0]) & (ZOY_CHdist_alt < ZOY_CHdist_cutoff[1]))) & \
                            (ZOY_NOdist > ZOY_NOdist_cutoff[0]) & (ZOY_NOdist < ZOY_NOdist_cutoff[1])] = 1
   
         if (len(listofZNcoordsH) > 0 and len(listofVcoordsO) > 0):
            ZNV_hbond_frame[(ZNV_OHdist > ZNV_OHdist_cutoff[0]) & (ZNV_OHdist < ZNV_OHdist_cutoff[1]) & \
                            (ZNV_NOdist > ZNV_NOdist_cutoff[0]) & (ZNV_NOdist < ZNV_NOdist_cutoff[1]) & \
                            (ZNV_HCdist > ZNV_HCdist_cutoff[0]) & (ZNV_HCdist < ZNV_HCdist_cutoff[1])] = 1
   
         if (len(listofZNcoordsH) > 0 and len(listofWcoordsO) > 0):
            ZNW_hbond_frame[(ZNW_OHdist > ZNW_OHdist_cutoff[0]) & (ZNW_OHdist < ZNW_OHdist_cutoff[1]) & \
                            (ZNW_ONdist > ZNW_ONdist_cutoff[0]) & (ZNW_ONdist < ZNW_ONdist_cutoff[1]) & \
                            (ZNW_HCdist > ZNW_HCdist_cutoff[0]) & (ZNW_HCdist < ZNW_HCdist_cutoff[1])] = 1
   
         if (len(listofZNcoordsH) > 0 and len(listofXcoordsO) > 0):
            ZNX_hbond_frame[(ZNX_OHdist > ZNX_OHdist_cutoff[0]) & (ZNX_OHdist < ZNX_OHdist_cutoff[1]) & \
                            (ZNX_ONdist > ZNX_ONdist_cutoff[0]) & (ZNX_ONdist < ZNX_ONdist_cutoff[1]) & \
                          (((ZNX_HCdist > ZNX_HCdist_cutoff[0]) & (ZNX_HCdist < ZNX_HCdist_cutoff[1])) | ((ZNX_HCdist > ZNX_HCdist_cutoff[0]) & (ZNX_HCdist < ZNX_HCdist_cutoff[1])))] = 1
   
         if (len(listofZNcoordsH) > 0 and len(listofYcoordsN) > 0):
            ZNY_hbond_frame[(ZNY_NHdist > ZNY_NHdist_cutoff[0]) & (ZNY_NHdist < ZNY_NHdist_cutoff[1]) & \
                            (ZNY_CHdist > ZNY_CHdist_cutoff[0]) & (ZNY_CHdist < ZNY_CHdist_cutoff[1]) & \
                          (((ZNY_NNdist > ZNY_NNdist_cutoff[0]) & (ZNY_NNdist < ZNY_NNdist_cutoff[1])) | ((ZNY_NNdist > ZNY_NNdist_cutoff[0]) & (ZNY_NNdist < ZNY_NNdist_cutoff[1])))] = 1
   
         Vs_frame = np.concatenate((UV_hbond_frame,ZOV_hbond_frame,ZNV_hbond_frame), axis=0)
         Ws_frame = np.concatenate((UW_hbond_frame,ZOW_hbond_frame,ZNW_hbond_frame), axis=0)
         Xs_frame = np.concatenate((UX_hbond_frame,ZOX_hbond_frame,ZNX_hbond_frame), axis=0)
         Ys_frame = np.concatenate((UY_hbond_frame,ZOY_hbond_frame,ZNY_hbond_frame), axis=0)

         current_frame = np.concatenate((Vs_frame,Ws_frame,Xs_frame,Ys_frame), axis=1)
   
         current_list = np.column_stack(np.nonzero(current_frame)).tolist()

         if current_list in whole_3darray:
            indexofhit = whole_3darray.index(current_list)
            numofhits[indexofhit] += 1
         else:
            whole_3darray.append(current_list)
            indexofhit = whole_3darray.index(current_list)
            numofhits[indexofhit] += 1

         os.remove(pdbnum)

print " "
print "There are " + str(len(whole_3darray)) + " unique frames of " + str(numframes - 1) + " total frames"
print " "

numofhits = numofhits[numofhits != 0]

whole_3dnumpyarray = np.array(whole_3darray)

np.save("arrayofinteractions",whole_3dnumpyarray)
np.save("numofhits",numofhits) 

np.savetxt("UVhbonds.csv", UV_hbond_array, delimiter=",", fmt="%i")
np.savetxt("UWhbonds.csv", UW_hbond_array, delimiter=",", fmt="%i")
np.savetxt("UXhbonds.csv", UX_hbond_array, delimiter=",", fmt="%i")
np.savetxt("UYhbonds.csv", UY_hbond_array, delimiter=",", fmt="%i")
np.savetxt("ZOVhbonds.csv", ZOV_hbond_array, delimiter=",", fmt="%i")
np.savetxt("ZOWhbonds.csv", ZOW_hbond_array, delimiter=",", fmt="%i")
np.savetxt("ZOXhbonds.csv", ZOX_hbond_array, delimiter=",", fmt="%i")
np.savetxt("ZOYhbonds.csv", ZOY_hbond_array, delimiter=",", fmt="%i")
np.savetxt("ZNVhbonds.csv", ZNV_hbond_array, delimiter=",", fmt="%i")
np.savetxt("ZNWhbonds.csv", ZNW_hbond_array, delimiter=",", fmt="%i")
np.savetxt("ZNXhbonds.csv", ZNX_hbond_array, delimiter=",", fmt="%i")
np.savetxt("ZNYhbonds.csv", ZNY_hbond_array, delimiter=",", fmt="%i")


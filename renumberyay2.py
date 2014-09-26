#!/usr/bin/python

import re

resnumber = 1
firstline = True

with open("yay2.pdb") as infile:
    with open("_yay2.pdb", 'w') as outfile:
        for line in infile:
            startATOM = re.search('^ATOM', line)
            if startATOM:
                if firstline:
                    lastoriginalnum = int(line[22:26])
                    firstline = False
                originalnum = int(line[22:26])
                if lastoriginalnum < originalnum:
                    resnumber += 1
                outfile.write("%s%4d%s" % (line[0:22],resnumber,line[26:80]))
                lastoriginalnum = originalnum
            else:
                outfile.write(line)


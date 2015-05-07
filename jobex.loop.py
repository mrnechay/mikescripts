#!/usr/bin/python

# call this script just as you would call jobex:
# jobex.loop.py -c 500 -ri -energy 5 -gcart 2 -gexp 2
# will call jobex with the same options, but will place it inside
# a loop which switches between cartesian and different kinds
# of internal coordinates in an attempt to converge geometry

import sys, grep

jobexOptions = " ".join(sys.argv[1:])



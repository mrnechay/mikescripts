#!/usr/bin/python

import sys

def boolgrep(*arg):
    found = False

    with open(arg[-1]) as file:
        for line in file:
            for ii in arg[0:-1]:
                if ii in line:
                    found = True

    return found

if __name__ == "__main__":
    print int(boolgrep(*sys.argv[1:]))

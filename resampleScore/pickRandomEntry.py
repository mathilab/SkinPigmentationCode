#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 10:37:29 2019
pickRandomEntry.py
@author: danju
"""

import sys
import getopt
import random

def main(argv):
    try:
        opts,args = getopt.getopt(argv, 'n:l:o:', ['num=','list=','out='])
    except getopt.GetoptError:
        sys.exit(2)
    for opt,arg in opts:
        if opt in ('-n','--num'):
            n = int(arg) # number of times to pick random entry in list
        elif opt in ('-l','--list'):
            path = arg # path to list of strings
        elif opt in ('-o','--out'):
            out = arg # export path

    # read list of strings
    f = open(path, 'r')
    l = f.readlines()

    # choose random entries in a list
    pick_list = [0] * n
    for x in range(n):
        i = random.randint(0,len(l)-1)
        pick_list[x] = l[i]

    # export picks as a list
    file = open(out, 'w')
    for line in pick_list:
        file.write(line)

if __name__ == '__main__':
    main(sys.argv[1:])

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 15:30:46 2018

@author: danju
"""

# Binary search of ordered list with duplicates
# Find index of lower and upper bound

import pandas as pd
import sys
import getopt
from IPython.utils.capture import capture_output

def main(argv):
    try:
        opts,args = getopt.getopt(argv, 's:l:h:', ['snp=','low=','high='])
    except getopt.GetoptError:
        sys.exit(2)
    for opt,arg in opts:
        if opt in ('-s','--snp'):
            snp_file = arg
        elif opt in ('-l','--low'):
            low = float(arg) # low bound
        elif opt in ('-h','--high'):
            high = float(arg)  # high bound

    df = pd.read_csv(snp_file, sep="\t", header=None)
    df.columns = ['chr','pos','fq']
    arr = df.loc[:,'fq']

    top = 0
    bot = len(df)
    # find index for lower bound
    if low <= 0:
        low_i = 1
    else:
        while top <= bot:
            mid = int( top + (bot - top)/2 )
            if float(arr[mid]) >= low:
                if float(arr[mid-1]) < low:
                    low_i = mid + 1
                    break
                else:
                    bot = mid - 1
            else:
                top = mid + 1

    top = 0
    bot = len(df)
    # find index for upper bound
    if high >= 1:
        high_i = len(df)
    else:
        while top <= bot:
            mid = int( top + (bot - top)/2 )
            if float(arr[mid]) <= high:
                if float(arr[mid+1]) > high:
                    high_i = mid + 1
                    break
                else:
                    top = mid + 1
            else:
                 if arr[mid-1]:
                    bot = mid - 1


    with capture_output() as c:
        print(low_i, high_i)
    c()

if __name__ == '__main__':
    main(sys.argv[1:])

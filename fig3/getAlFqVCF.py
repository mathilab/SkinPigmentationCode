#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 15:33:27 2019

Calculate allele frequency of alternate alleles from VCF and output CSV

@author: danju
"""

import csv
import gdc
import sys
import getopt

def tally(count,geno,d,total):
    if d == True:
        count = count + int(geno[0]) + int(geno[2])
        total = total + 2
    else:
        count = count + int(geno[0])
        total = total + 1
    return count, total

def main(argv):
    try:
        opts,args = getopt.getopt(argv, "o:v:d:", ["output=","vcf=","diploid="])
    except getopt.GetoptError:
        sys.exit(2)
    for opt,arg in opts:
        if opt in ("-o","--output"):
            output_filename = arg
        elif opt in ("-v","--vcf"):
            vcf_filename = arg
        elif opt in ("-d","--diploid"): # samples diploid? (True/False)
            if (arg in ['1', 1, 'T', 't', 'True', 'true']):
                d_flag = True
            else:
                d_flag = False

    n_head_lines = 0
    for vcfline in gdc.open2(vcf_filename):
        n_head_lines += 1
        line = vcfline.decode()
        if line[1] != '#':
            break
    nlines = sum(1 for line in gdc.open2(vcf_filename)) - n_head_lines

    vcf_file = gdc.open2(vcf_filename)
    vcfline=next(vcf_file)
    pound = vcfline.decode()[1]
    while pound=='#':
        vcfline=next(vcf_file)
        pound = vcfline.decode()[1]
    line = vcfline.decode()

    output_head = ['ID','chr','pos','Alt','N','Fq']

    with open(output_filename,'w') as mycsv:
        datawriter = csv.writer(mycsv)
        datawriter.writerow(output_head)
#, quoting = csv.QUOTE_NONE
        for i in range(0,nlines):
            vcfline = next(vcf_file)
            line = vcfline.decode()
            line = line.split()
            chr, pos = line[0], line[1]
            chrpos = line[0] + '_' + line[1]

            alt_count = 0
            total = 0
            for j in range(9,len(line)):
                if line[j][0] == '.':
                    continue
                else:
                    alt_count, total = tally(alt_count,line[j],d_flag,total)

            if total == 0:
                continue
            else:
                alt_fq = alt_count / total
                row = [chrpos,chr,pos,alt_count,total,alt_fq]
                datawriter.writerow(row)

    mycsv.close()

if __name__ == "__main__":
    main(sys.argv[1:])
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 09:13:08 2018

@author: danju
"""

#Returns a matrix from a VCF of whether sample has alternative SNP
 #Handles biallelic SNPs only

import csv
import sys
import getopt
import gzip


def main(argv):
    try:
        opts,args = getopt.getopt(argv, 'o:v:a:', ['output=','vcf=','apulldown='])
    except getopt.GetoptError:
        sys.exit(2)
    for opt,arg in opts:
        if opt in ('-o','--output'):
            output_filename = arg
        elif opt in ('-v','--vcf'):
            vcf_filename = arg
        elif opt in ('-a','--apulldown'): # VCF format from apulldown?
            apull_format = arg

    n_head_lines = 0
    for vcfline in gzip.open(vcf_filename):
        n_head_lines += 1
        line = vcfline.decode()
        if line[1] != '#':
            break

    nlines = sum(1 for line in gzip.open(vcf_filename)) - n_head_lines

    samples = ['SNP']
    vcf_file = gzip.open(vcf_filename)
    vcfline=next(vcf_file)
    pound = vcfline.decode()[1]
    while pound=='#':
        vcfline=next(vcf_file)
        pound = vcfline.decode()[1]
    line = vcfline.decode()
    header = line.split('\t')
    for i in range(9,len(header)):
        sample_name = header[i]
        if '\n' in sample_name: #handle last sample in header
            sample_name = sample_name[0:-1]
        samples.append(sample_name)

#    line=next(vcf_file)
#    line=line.decode()
    with open(output_filename,'w') as mycsv:
        datawriter = csv.writer(mycsv)
        datawriter.writerow(samples) #header line

        for i in range(0,nlines):
            chr_pos = []
            vcfline=next(vcf_file)
            line=vcfline.decode()
            line = line.split()
            cur_chrpos = line[0] + '_' + line[1]
            chr_pos.append(cur_chrpos)
            geno_list = line[9:]

            # handle apulldown vcf format
            if apull_format=='T':
               alt_str = ''
               for j in geno_list:
                   alt = j[0]
                   alt_str += alt
               out_line = chr_pos + list(alt_str)
            else:
                # captureShotgun format
                geno_str = ''.join(geno_list)
                geno = list(geno_str[0::3])
                out_line = chr_pos + geno

            datawriter.writerow(out_line)

    mycsv.close()

if __name__ == '__main__':
    main(sys.argv[1:])

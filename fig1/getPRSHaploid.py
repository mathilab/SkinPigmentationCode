#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 12:33:45 2018

@author: danju
"""
## Calculate pigmentation polygenic risk score from VCFs and a 'master' file
#` of variants with direction of effect alleles and betas
## Update 9/12: weighted PRS is now calculated as a fraction of weighted value

import csv
import random
import sys
import getopt

from countDarkAlleles import getNumDenom
from calcPRS import getPRS

def main(argv):
    try:
        opts,args = getopt.getopt(argv, 'o:v:m:b:e:', ['output=','vcflist=',
                                                       'master=','use_beta=',
                                                       'eff_allele='])
    except getopt.GetoptError:
        sys.exit(2)
    for opt,arg in opts:
        if opt in ('-o','--output'): # output file name
            output_filename = arg
        elif opt in ('-v','--vcflist'): # each row a VCF file name and haploid
                                        # status (1/0). If 0, diploid samples
                                        # are separated into two haploids.
            vcf_list_filename = arg
        elif opt in ('-m','--master'): # CSV file name with SNPs in proper
                                       # format
            master_filename = arg
        elif opt in ('-b','--use_beta'): # use effect size in PRS calculation?
            use_beta = arg
            if (use_beta in ['1', 1, 'T', 't', 'True', 'true']):
                    use_beta = True
            else:
                use_beta = False
        elif opt in ('-e','--eff_allele'):
            eff_col = int(arg)

    with open(master_filename,encoding='utf8') as mastercsv:
        masterreader = csv.reader(mastercsv)
        snp_dict = {}
        for snp in masterreader:
            chrom_pos = snp[0] + '_' + snp[1]
            dark, light = snp[3], snp[4]
            eff, beta = snp[eff_col], snp[11]
            temp_entry_snp = {chrom_pos:[dark,light,eff,beta]}
            snp_dict.update(temp_entry_snp)

    # make dictionary of vcf files and whether samples have haploid calls (0/1)
    vcf_file_list = open(vcf_list_filename,'r')
    vcf_dict = {}
    for vcf_file in vcf_file_list:
        vcf_file = vcf_file.split('\t')
        vcfname, haploid_status = vcf_file[0], vcf_file[1]
        temp_entry_vcf = {vcfname:haploid_status}
        vcf_dict.update(temp_entry_vcf)

    sample_score_dict = {}
    for vcf_filename in vcf_dict.keys(): # loop through VCFs
        cur_hap_status = int(vcf_dict[vcf_filename])
        vcf_file = open(vcf_filename,'r')
        for vcfline in vcf_file: # loop through SNPs
            if vcfline[:2] == '##':
                continue

            if vcfline[0]=='#' and vcfline[1] != '#':
                header = vcfline.split('\t')
                # make dictionary for samples:[dark,total,prs,weight]
                # note: weight only used when using gwas eff size and is
                # the number of SNPs present in sample
                for i in range(9,len(header)):
                    sample_name = header[i]
                    if '\n' in sample_name: # handle last sample in header
                        sample_name = sample_name[0:-1]
                    if cur_hap_status == 1:
                        temp_entry_score = {sample_name:[0,0,0,0]}
                        sample_score_dict.update(temp_entry_score)
                    elif cur_hap_status == 0:
                        name1 = sample_name + '_chr1'
                        name2 = sample_name + '_chr2'
                        temp_entry_score1 = {name1:[0,0,0,0]}
                        temp_entry_score2 = {name2:[0,0,0,0]}
                        sample_score_dict.update(temp_entry_score1)
                        sample_score_dict.update(temp_entry_score2)

            if vcfline[0] != '#':
                vcfline = vcfline.split('\t')
                cur_chrpos = vcfline[0] + '_' + vcfline[1]
                cur_ref, cur_alt = vcfline[3], vcfline[4]
                cur_dark = snp_dict[cur_chrpos][0]
                cur_eff = snp_dict[cur_chrpos][2]
                cur_beta = float(snp_dict[cur_chrpos][3])
#                list_miss = [] # list of samples with missing genotype for SNP
#                num_geno_sample = 0 # counter for samples with genotype data
#                cur_total_beta = 0 # running sum of betas

                for i in range(9,len(header)): # loop through samples
                    cur_sample = header[i]
                    if '\n' in cur_sample: # handle last sample in curent line
                        cur_sample = cur_sample[0:-1]
                    genotype = vcfline[i]
                    if genotype[0] == '.' or genotype[2] == '.': # no genotype
#                        list_miss.append(cur_sample)
                        continue
                    # randomize which allele so no bias for diploid samples
                    ran_int = random.randint(0,1)
                    if ran_int == 1:
                        allele1, allele2 = int(genotype[0]), int(genotype[2])
                    else:
                        allele1, allele2 = int(genotype[2]), int(genotype[0])
                    # compute PRS
                    if cur_hap_status == 1: # haploid
                        if use_beta:
                            old_prs = sample_score_dict[cur_sample][0]
                            try:
                                cur_prs1 = getPRS(old_prs,allele1,cur_alt,cur_ref,
                                                 cur_eff,cur_beta)
                                cur_prs2 = getPRS(old_prs,allele2,cur_alt,cur_ref,
                                                 cur_eff,cur_beta)
                                cur_prs = (cur_prs1 + cur_prs2) / 2
                            except UnboundLocalError:
                                print('Effect allele not matching genotype.\
                                      Check master SNP file.' )
                                print(cur_chrpos, cur_sample)
                                break
                            sample_score_dict[cur_sample][0] = cur_prs
                            sample_score_dict[cur_sample][1] += abs(cur_beta)
                            sample_score_dict[cur_sample][3] += 1
#                            # build average for a SNP (now deprecated)
#                            cur_sample_beta = cur_prs - old_prs
#                            cur_total_beta = cur_total_beta + cur_sample_beta
#                            num_geno_sample = num_geno_sample + 1

                        elif not use_beta:
                            old_num = sample_score_dict[cur_sample][0]
                            old_denom = sample_score_dict[cur_sample][1]
                            try:
                                cur_num, cur_denom = getNumDenom(old_num,old_denom,
                                                                   allele1,cur_alt,
                                                                   cur_dark,cur_ref)
                            except UnboundLocalError:
                                print('Effect allele not matching genotype.\
                                      Check master SNP file.' )
                                print(cur_chrpos, cur_sample)
                                break
                            sample_score_dict[cur_sample][0] = cur_num
                            sample_score_dict[cur_sample][1] = cur_denom

                    elif cur_hap_status == 0: # diploid
                        cur_sam1 = cur_sample + '_chr1'
                        cur_sam2 = cur_sample + '_chr2'

                        if use_beta:
                            old_prs1 = sample_score_dict[cur_sam1][0]
                            old_prs2 = sample_score_dict[cur_sam2][0]
                            cur_prs1 = getPRS(old_prs1,allele1,cur_alt,cur_ref,
                                              cur_eff,cur_beta)
                            cur_prs2 = getPRS(old_prs2,allele2,cur_alt,cur_ref,
                                              cur_eff,cur_beta)
                            sample_score_dict[cur_sam1][0] = cur_prs1
                            sample_score_dict[cur_sam2][0] = cur_prs2
                            sample_score_dict[cur_sam1][1] += abs(cur_beta)
                            sample_score_dict[cur_sam2][1] += abs(cur_beta)
#                            cur_sample_beta1 = cur_prs1 - old_prs1
#                            cur_sample_beta2 = cur_prs2 - oEPRECATED)
#                            cur_total_beta = (cur_total_betald_prs2
#                            # build average for a SNP (NOW D + cur_sample_beta1
#                                             + cur_sample_beta2)
#                            num_geno_sample = num_geno_sample + 2

                        elif not use_beta:
                            old_num1 = sample_score_dict[cur_sam1][0]
                            old_denom1 = sample_score_dict[cur_sam1][1]
                            old_num2 = sample_score_dict[cur_sam2][0]
                            old_denom2 = sample_score_dict[cur_sam2][1]
                            cur_num1, cur_denom1 = getNumDenom(old_num1,old_denom1,
                                                                 allele1,cur_alt,
                                                                 cur_dark,cur_ref)
                            cur_num2, cur_denom2 = getNumDenom(old_num2,old_denom2,
                                                                 allele2,cur_alt,
                                                                 cur_dark,cur_ref)
                            sample_score_dict[cur_sam1][0] = cur_num1
                            sample_score_dict[cur_sam1][1] = cur_denom1
                            sample_score_dict[cur_sam2][0] = cur_num2
                            sample_score_dict[cur_sam2][1] = cur_denom2

#                if use_beta: (DEPRECATED)
#                    if num_geno_sample == 0:
#                        break
#                    snp_avg = cur_total_beta / num_geno_sample
#                    for name in list_miss: # impute score for SNP
#                        sample_score_dict[name][2] = (sample_score_dict[name][2]
#                                                     + snp_avg)

        vcf_file.close()

    # write to output file
    with open(output_filename,'w') as mycsv:
        datawriter = csv.writer(mycsv)
        if use_beta:
            datawriter.writerow(['ID','Score','NumDark','TotalNumSNPs','Weight'])
        elif not use_beta:
            datawriter.writerow(['ID','Score','NumDark','TotalNumSNPs'])
        sample_list = list(sample_score_dict.keys())
        for key in sample_list:
            if use_beta:
                num = sample_score_dict[key][0]
                denom = sample_score_dict[key][1] # highest possible score
                if denom == 0:
                    continue
                score = num / denom
                weight = sample_score_dict[key][3]
                datawriter.writerow([key,score,num,denom,weight])
            elif not use_beta:
                denom = sample_score_dict[key][1]
                if denom == 0:
                    continue
                num = sample_score_dict[key][0]
                score = num / denom
                datawriter.writerow([key,score,num,denom])

    mycsv.close()

if __name__ == '__main__':
    main(sys.argv[1:])

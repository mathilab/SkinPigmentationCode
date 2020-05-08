#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 13:37:30 2018

@author: danju
"""
# Calculate polygenic risk score based on proportion of dark alleles

def getNumDenom(old_num,old_denom,allele,curr_alt,curr_dark,curr_ref):
    if len(curr_alt) > 1: 
        curr_alt1, curr_alt2 = curr_alt[0], curr_alt[2]
    else:
        curr_alt1, curr_alt2 = curr_alt[0], 'NA' 
    
    if curr_alt1 == curr_dark:
        if allele == 1:
            curr_num = old_num + 1
        else:
            curr_num = old_num
            
    elif curr_alt2 == curr_dark:
        if allele == 2:
            curr_num = old_num + 1
        else: 
            curr_num = old_num
            
    elif curr_ref == curr_dark:
        if allele == 0:
            curr_num = old_num + 1
        else: 
            curr_num = old_num    
    curr_denom = old_denom + 1
    
    return curr_num, curr_denom
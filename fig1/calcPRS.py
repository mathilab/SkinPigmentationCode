#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 16:51:36 2018

@author: danju
"""
# Calculate polygenic risk score for skin color using effect size
# Polarizes effect allele so beta always in direction of dark (+)

def getPRS(old_prs,allele,curr_alt,curr_ref,curr_eff,beta):
    # handle multiallelic sites
    if len(curr_alt) > 1:
        curr_alt1, curr_alt2 = curr_alt[0], curr_alt[2]
    else:
        curr_alt1, curr_alt2 = curr_alt[0], 'NA'

    if curr_alt1 == curr_eff:
        if beta > 0:
            if allele == 1:
                curr_prs = old_prs + beta
            else:
                curr_prs = old_prs
        else:
            if allele == 0:
                curr_prs = old_prs + abs(beta)
            else:
                curr_prs = old_prs

    elif curr_alt2 == curr_eff:
        if beta > 0:
            if allele == 2:
                curr_prs = old_prs + beta
            else:
                curr_prs = old_prs
        else:
            if allele == 0:
                curr_prs = old_prs + abs(beta)
            else:
                curr_prs = old_prs

    elif curr_ref == curr_eff:
        if beta > 0:
            if allele == 0:
                curr_prs = old_prs + beta
            else:
                curr_prs = old_prs
        else:
            if allele == 1 or allele == 2:
                curr_prs = old_prs + abs(beta)
            else:
                curr_prs = old_prs

    return curr_prs
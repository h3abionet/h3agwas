#!/usr/bin/env python

from __future__ import print_function


import sys
import re

f = open(sys.argv[1])
rem_name=sys.argv[2]
rem_inds = rem_snps = rem_maf = rem_hwe =  0

autosome=False
nonautosome = 0
for line in f:
    if "Options in effect" in line: 
        autosome=False
        continue
    if "--autosome" in line:
        autosome=True
        continue
    m =re.search("(\w+) out of (\w+) variants loaded from .bim file.", line)
    if m and autosome:
        nonautosome = int(m.group(2))-int(m.group(1))
        continue
    m=re.search("(\w+) (\w+) removed due to ([-A-z0-9]+ \w+)",line)
    if m:
        #print(m.group(1),m.group(2),m.group(3))
        n = int(m.group(1))
        if m.group(2)=="variants":
            if  m.group(3)=="missing genotype":
                rem_snps = rem_snps+n
            elif m.group(3)=="Hardy-Weinberg exact":
                rem_hwe = rem_hwe+n
            else:
                rem_maf = rem_maf + n
        else:
            rem_inds = rem_inds + n

text = """
*-noindent
Using this approach, 
*-begin{itemize}
*-item %d SNPs that are non-autosomal were removed;
*-item %d SNPs were removed due missing genotype threshold constraints;
*-item %d individuals were removed due to missing genotype constraints (the list of missing individuals, if any, can be found in the file *-url{%s});
*-item %d SNPs were removed as the MAF was too low.
*-item %d SNPs were removed as they were out of the specified Hardy-Weinberg equilibrium. 
*-end{itemize}
""" %(nonautosome,rem_snps,rem_inds,rem_name,rem_maf,rem_hwe)

print(text.replace("*-",chr(92)).replace("##",chr(36)))

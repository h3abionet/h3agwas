#!/usr/bin/env python3

# we allow  the reference file to be in multiple format 
# 1. Simple: 
#  no header, two columns first column is the SNP name, second is the reference allele
# 2. complex
#  as many columns as you like -- the first column is a header, one column name is SNP one column is Base 
#  these are the only columns that matter (to this workflow)


import sys
import os

inp = sys.argv[1]
out = sys.argv[2]

fields = open(inp).readline().strip().split()

if "empty" in inp: 
    print ("")
    sys.exit(0)

if "SNP" not in fields and "Base" not in fields:
    if len(fields)==2:
        parm = "2 1"
    else:
        sys.exit("Illegal file format in <%s> (too many columns)"%inp)
elif "SNP" in fields and "Base" in fields:
    snp_col = fields.index("SNP")+1
    all_col = fields.index("Base")+1
    parm = "%d %d 1 "%(all_col,snp_col)
else:
   sys.exit("Illegal file format in <%s> (header not legal)"%inp)

print(parm)



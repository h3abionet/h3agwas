#!/usr/bin/env python3
from __future__ import print_function
import sys

# Called as a template from nextflow
if len(sys.argv) == 1:
    sys.argv=["dups.py","$inpfname","$outfname","$remove_on_bp"]


def removeOnBP(fname, out):
    f = open(fname)
    line = f.readline().strip().split()
    old_chrom = line[0]
    old_snp   = line[1]
    old_bp    = line[2]
    for line in f:
        data = f.readline().strip().split()
        chrom = data[0]
        snp   = data[1]
        bp    = data[2]
        if (chrom, bp) == (old_chrom, old_bp):
            out.write("%s\\n"%snp)
        elif (chrom, bp) > (old_chrom, old_bp):
            (old_chrom,old_snp,old_bp) = (chrom,snp,bp)
        else:
            print(""" 

              The BIM file <%s> is not sorted see <%s> and <%s>"
            """ % (fname, old_snp, snp))



f=open(sys.argv[1])
if not f:
    sys.exit("File <%s> not opened"%sys.argv[1])
s_name = set()
out=open(sys.argv[2],"w")
for line in f:
    data=line.strip().split()
    snp_name = data[1]
    if snp_name in s_name:
        out.write( "%s\\n"%snp)
    else:
        s_name.add(data[1])
f.close()

if sys.argv[3] in ["1",1,True,"True","true"]:
    removeOnBP(sys.argv[1],out)

out.close()


#!/usr/bin/env python3
from __future__ import print_function
import sys
import os

# Called as a template from nextflow
if len(sys.argv) == 1:
    sys.argv=["dups.py","$inpfname","$outfname","$remove_on_bp"]


def getChrom(chrom):
    try:
        result = int(chrom)
    except ValueError:
        result = chrom
    return result

def removeOnBP(fname, out):
    f = open(fname)
    line = f.readline().strip().split()
    old_chrom = getChrom(line[0])
    old_snp   = line[1]
    old_bp    = int(line[3])
    for line in f:
        data = line.strip().split()
        try:
            chrom = getChrom(data[0])
        except IndexError:
            print("****",line)
            print(data)
            print(old_snp)
            sys.exit(11)
        snp   = data[1]
        bp    = int(data[3])
        if (chrom, bp) == (old_chrom, old_bp):
            out.write("%s\\n"%snp)
            (old_chrom,old_snp,old_bp) = (chrom,snp,bp)
        elif (chrom, bp) > (old_chrom, old_bp):
            (old_chrom,old_snp,old_bp) = (chrom,snp,bp)
        else:
            print(old_chrom, old_bp, old_snp)
            print(chrom,bp,snp)
            print(""" 

              The BIM file <%s> is not sorted see <%s> and <%s>"
            """ % (fname, old_snp, snp))
            sys.exit(10)


os.system("hostname > hostname")
f=open(sys.argv[1])
if not f:
    sys.exit("File <%s> not opened"%sys.argv[1])
s_name = set()
out=open(sys.argv[2],"w")
for line in f:
    data=line.strip().split()
    snp_name = data[1]
    if snp_name in s_name:
        out.write( "%s\\n"%snp_name)
    else:
        s_name.add(data[1])
f.close()

if sys.argv[3] in ["1",1,True,"True","true"]:
    removeOnBP(sys.argv[1],out)

out.close()


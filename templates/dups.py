#!/usr/bin/env python3
from __future__ import print_function
import sys

# Called as a template from nextflow
if len(sys.argv) == 1:
    sys.argv=["dups.py","$inpfname","$outfname"]

f=open(sys.argv[1])
if not f:
    sys.exit("File <%s> not opened"%sys.argv[1])
s=set()
out=open(sys.argv[2],"w")
for line in f:
    data=line.strip().split()
    snp=data[1]
    if snp in s:
        out.write( "%s\\n"%snp)
    s.add(data[1])
f.close()
out.close()


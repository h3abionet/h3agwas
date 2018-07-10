#!/usr/bin/env python3

import sys

""" 
Take a relatedness file from GEMMA and tranform into version for FASTLMM. 
We need a header line with the IDs of the individuals (FID IID separated by space) separated by TAB.
Then each line has as first column the FID IID pair followed by the distances to all the other individuals.
"""


f        = open(sys.argv[1])
rel_base = open(sys.argv[2])
outf     = open(sys.argv[3],"w")

my_ids  = []
for line in f:
    data = line.split()
    new_id = "%s %s"%(data[0],data[1])
    my_ids.append(new_id)
outf.write("var\t"+("\t".join(my_ids))+"\n")
f.close()


for curr in my_ids:
    rel = rel_base.readline()
    rel.replace(" ","\t")
    outf.write(curr+"\t"+rel+"\n")
rel_base.close()

outf.close()
    



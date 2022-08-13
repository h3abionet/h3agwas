#!/usr/bin/env python3
import hashlib
import argparse

from os.path import basename

def parseArguments():
    parser = argparse.ArgumentParser(description='plot pca using python')
    parser.add_argument('--bfile',type=str,required=True,help="File with phenotype and covariate data")
    parser.add_argument('--out',type=str,required=True,help="case control column")
    args = parser.parse_args()
    return args


# Taken from http://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file

def md5sum(filename, blocksize=65536):
    hash = hashlib.md5()
    with open(filename, "rb") as f:
        for block in iter(lambda: f.read(blocksize), b""):
            hash.update(block)
    return hash.hexdigest()

args = parseArguments()


md5s=[]
plinks=[ args.bfile+".bed",args.bfile+".bim",args.bfile+".fam"]

for x in plinks:
    md5s.append((basename(x),md5sum(x)))




f=open(args.out,"w")
for (fn,md) in md5s:
    f.write("%s\t%s\n"%(fn,md))
f.close()


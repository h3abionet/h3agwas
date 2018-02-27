#!/usr/bin/env python3


import pandas as pd
import argparse
import numpy as np
import sys

# Preferably a genuine strand file is given -- but if a dummy strandfile is give then
# we use the manifest file -- it is much much slower


def parseArguments():
    parser = argparse.ArgumentParser(description='fill in missing bim values')
    parser.add_argument('strand',type=str)
    parser.add_argument('manifest',type=str)
    parser.add_argument('bim', type=str)
    parser.add_argument('output', type=str)
    args = parser.parse_args()
    return args


TAB=chr(9)
EOL=chr(10)

def rc(x):
    if x=='A': return 'T'
    elif x=='C': return 'G'
    elif x=='G': return 'C'
    elif x=='T': return 'A'
    return x

def strandProc(bdata):
     alleles = mf.loc[bdata[1]].values
     return alleles

def manifestProc(bdata):
     data = mf.loc[bdata[1]].values
     alleles=data[1][1:-1].split("/")
     if data[0] == "BOT":
         alleles = list(map(rc,alleles))
     return alleles

args = parseArguments()


if args.strand == "emptyZ0strn.txt": # we have a dummy strand file and  the manifest file is given
    mf = pd.read_csv(args.manifest,delimiter=",",skiprows=7,\
                     usecols=["SNP","IlmnStrand","Name"])
    name='Name'
    extractValues = manifestProc
else:
    mf = pd.read_csv(args.strand,delim_whitespace=True,\
                     usecols=["SNP_Name","Top_AlleleA", "Top_AlleleB"],comment="#")
    name = "SNP_Name"
    extractValues = strandProc
mf.set_index(name,inplace=True)


bf = open(args.bim)

out=open(args.output,"w")
for line in bf:
    bdata = line.rstrip().split()
    if bdata[4] == '0':
        alleles = extractValues(bdata)
        if bdata[5]   == alleles[0]:
            bdata[4]=alleles[1]
        elif bdata[5] == alleles[1]:
            bdata[4]=alleles[0]
        elif bdata[5]=='0':
            bdata[4]=alleles[0]
            bdata[5]=alleles[1]
        else:
            print("No good minor allele {}  {}".format(line.strip(),str(alleles)))
    out.write(TAB.join(bdata)+EOL)
out.close()


                                
         




#!/usr/bin/env python3

# (c) University of the Witwatersand, Johannesburg on behalf of the H3ABioinformatics Network Consortium
# 2016-2018
# Licensed under the Creative Commons Attribution 4.0 International Licence. 
# See the "LICENSE" file for details

import pandas as pd
import argparse
import numpy as np
import sys

# Preferably a genuine strand file is given -- but if a dummy strandfile is give then
# we use the manifest file -- it is much much slower


def parseArguments():
    parser = argparse.ArgumentParser(description='fill in missing bim values')
    parser.add_argument('align',type=str)
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


if args.align == "db2ref":
    all1_col, all2_col  =   "Forward_Allele1", "Forward_Allele2"
else:
    all1_col, all2_col  =   "Top_AlleleA", "Top_AlleleB"

if args.strand == "emptyZ0strn.txt": # we have a dummy strand file and  the manifest file is given
    mf = pd.read_csv(args.manifest,delimiter=",",skiprows=7,\
                     usecols=["SNP","IlmnStrand","Name"])
    name='Name'
    extractValues = manifestProc
else:
    print("am here",all1_col)
    mf = pd.read_csv(args.strand,delim_whitespace=True,\
                     usecols=["SNP_Name",all1_col,all2_col],comment="#")
    name = "SNP_Name"
    extractValues = strandProc
mf.set_index(name,inplace=True)


bf = open(args.bim)

out=open(args.output,"w")
for line in bf:
    bdata = line.rstrip().split()
    if "rs99365" in bdata[1]: 
        print(bdata)
        alleles = extractValues(bdata)
        print(alleles)
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


                                
         




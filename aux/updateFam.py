#!/usr/bin/env python3

import argparse
import pandas as pd
import sys
import re

def parseArguments():
    parser=argparse.ArgumentParser()
    parser.add_argument('samplesheet', type=str, metavar='samplesheet'),
    parser.add_argument('updatesheet', type=str, metavar='updatesheet'),
    parser.add_argument('oldfam', type=str, metavar='oldfam'),
    parser.add_argument('newfam', type=str, metavar='newfam'),    
    args = parser.parse_args()
    if args.newfam == args.oldfam:
        sys.exit("New and old fam cannot be the same")
    return args


PIDS = ['FID','IID']
 
def getFam(famf):
   fam =  pd.read_csv(famf,delim_whitespace=True,header=None,\
                             names=["FID","IID","FAT","MAT","SEX","PHE"])
   old = ['OFID','OIID']
   fam['OFID']=fam['FID']
   fam['OIID']=fam['IID']
   fam.set_index(old, inplace=True)
   return fam
      
def getSmplLbl(data):
    m=re.search(".*_(\w+)",data["Institute Sample Label"])
    if not m: sys.exit("Failed parsing "+data)
    return m.group(1)



args = parseArguments()
fam = getFam(args.oldfam)

orig = pd.read_excel(args.samplesheet)
orig.set_index(["Institute Plate Label","Well"],inplace=True)
update = pd.read_excel(args.updatesheet,skiprows=3)
count=0
for i, row in update.iterrows():
    if "IlluminaControl" in row["Institute Sample Label"]: continue
    new_lbl=getSmplLbl(row)
    pos      = tuple(row[['Institute Plate Label','Well']].values)
    oldid    = getSmplLbl(orig.loc[pos])
    try:
        if oldid != new_lbl:
            print("%s --> %s"%(oldid,new_lbl),len(fam))
        fam.loc[(oldid,oldid)]['FID']   = new_lbl
        fam.loc[(oldid,oldid)]['IID']   = new_lbl
    except KeyError:
        print("No data for "+new_lbl+"!!!!")
        continue
    if oldid != new_lbl:
        count=count+1

    
fam.to_csv(args.newfam,sep="\t",header=None,index=False),

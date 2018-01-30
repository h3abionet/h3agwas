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
    parser.add_argument('output', type=str, metavar='output'),    
    args = parser.parse_args()
    if "%s.fam"%args.output == args.oldfam:
        sys.exit("New and old fam cannot be the same")
    return args

EOL=chr(10)
PIDS = ['FID','IID']
 
def getFam(famf):
   fam =  pd.read_csv(famf,delim_whitespace=True,header=None,\
                             names=["FID","IID","FAT","MAT","SEX","PHE"])
   old = ['OFID','OIID']
   fam['OFID']=fam['FID']
   fam['OIID']=fam['IID']
   fam.set_index(old, inplace=True)
   oldfam = fam.copy(deep=True)
   return oldfam,fam
      
def getSmplLbl(data):
    m=re.search(".*_(\w+)",data["Institute Sample Label"])
    if not m: sys.exit("Failed parsing "+data)
    return m.group(1)



args = parseArguments()
oldfam,fam = getFam(args.oldfam)

orig = pd.read_excel(args.samplesheet)
orig.set_index(["Institute Plate Label","Well"],inplace=True)
orig['PID']=orig.apply(getSmplLbl,axis=1)
update = pd.read_excel(args.updatesheet,skiprows=3)
count=0
g=open("%s.errs"%args.output,"w")
h=open("%s.switch"%args.output,"w")

for i, row in update.iterrows():
    if "IlluminaControl" in row["Institute Sample Label"]: continue
    new_lbl=getSmplLbl(row)
    pos      = tuple(row[['Institute Plate Label','Well']].values)
    oldid    = getSmplLbl(orig.loc[pos])
    if oldid != new_lbl:
        h.write("%s  -> %s \n"%(oldid,new_lbl))
        count=count+1
        if fam.index.contains((new_lbl,new_lbl)):
            fam.loc[(oldid,oldid),["FID","IID"]]   = new_lbl
            fam.loc[(oldid,oldid),["FAT","MAT","SEX","PHE"]]   = oldfam.loc[(new_lbl,new_lbl),["FAT","MAT","SEX", "PHE"]]
        else:
            g.write(new_lbl+EOL)
            continue
    if "unknown" in row["Institute Sample Label"]: 
        print("Pos=<%s>; Old =<%s>; New=<%s>"%(pos,oldid,new_lbl))
        print(fam.loc[(oldid,oldid)]['IID'])


g.close()    
h.close()
all_ids = fam[PIDS].to_records(index=False).tolist()
uniq =set([])
orig.set_index('PID',inplace=True)
for x in all_ids:
    if x in uniq:
        print("Duplicated element ",x[0],orig.loc[x[0],"Institute Sample Label"])
    else:
        uniq.add(x)

fam.to_csv("%s.fam"%args.output,sep="\t",header=None,index=False),

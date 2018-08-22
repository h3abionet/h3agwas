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


def sex_code(x):
    if x in  ["Male","M"]:
        return "1"
    elif x in  ["Female","F"]:
        return "2"
    else:
        return "0"

 
def getFam(famf):
   fam =  pd.read_csv(famf,delim_whitespace=True,header=None,\
                             names=["FID","IID","FAT","MAT","SEX","PHE"],\
                             dtype={"FAT":str,"MAT":str,"SEX":str, "PHE":str,"FID":str,"IID":str})
   old = ['OFID','OIID']
   fam['OFID']=fam['FID']
   fam['OIID']=fam['IID']
   fam.set_index(old, inplace=True)
   oldfam = fam.copy(deep=True)
   return oldfam,fam
      
def getSmplLbl(data):
    m=re.search(".*_(\w+)",data["Institute Sample Label"])
    if not m: sys.exit("Failed parsing "+data)
    return m.group(1).replace(" ","")



args = parseArguments()
oldfam,fam = getFam(args.oldfam)

orig = pd.read_excel(args.samplesheet,dtype={"Institute Sample Label":str})
orig.set_index(["Institute Plate Label","Well"],inplace=True)
orig['PID']=orig.apply(getSmplLbl,axis=1)


def performUpdate(g,h, sheet):
    update = pd.read_excel(sheet,skiprows=3)
    m = re.search(r"batch(\d+)",sheet)
    if m:
        bnumber = m.group(1)
    else:
        sys.exit("Can't get batch number in file name")
    count=0
    for i, row in update.iterrows():
        if "IlluminaControl" in row["Institute Sample Label"]: continue
        new_lbl=getSmplLbl(row)    # label in the update sheet
        pos      = tuple(row[['Institute Plate Label','Well']].values) # where it was plate/well
        oldid    = getSmplLbl(orig.loc[pos])  # What that plate/well was labelled as originally
        if oldid != new_lbl:  # If there's a change
            h.write("%s  -> %s \n"%(oldid,new_lbl))
            count=count+1
            # Was this in the original FAM file
            if oldfam.index.contains((oldid,oldid)):
                fam.loc[(oldid,oldid),["FID","IID"]]   = new_lbl
                if oldfam.index.contains((new_lbl,new_lbl)):
                    if new_lbl  == "BFW0K": print (1,oldid,new_lbl)
                    for col in ["FAT","MAT","SEX","PHE"]:
                        value = oldfam.loc[(new_lbl,new_lbl),col]
                        if type(value) != str:
                            value = value.values[0]
                        fam.loc[(oldid,oldid),col]   = value
                else:
                    fam.loc[(oldid,oldid),"SEX"]=sex_code(row["Sex"])
                    fam.loc[(oldid,oldid),"MAT"]=0
                    fam.loc[(oldid,oldid),"FAT"]=0
                    fam.loc[(oldid,oldid),"PHE"]=bnumber
            else:
                g.write(oldid+" unknown in original fam file "+EOL)
                continue
        if "unknown" in row["Institute Sample Label"]: 
           g.write("unknown label at %s (%s %s  %s)"%(pos,oldid,new_lbl,fam.loc[(oldid,oldid)]['IID']))


g=open("%s.errs"%args.output,"w")
h=open("%s.switch"%args.output,"w")
for sheet in args.updatesheet.split(","):
    performUpdate(g,h,sheet)
g.close()    
h.close()
all_ids = fam[PIDS].to_records(index=False).tolist()
uniq =set([])
orig.set_index('PID',inplace=True)
for x in all_ids:
    if x in uniq:
        print("Duplicated element ",x[0])
        pass
    else:
        uniq.add(x)

fam.to_csv("%s.fam"%args.output,sep="\t",header=None,index=False),

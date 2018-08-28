#!/usr/bin/env python3

"""
Input:
-- FAM file :   the fam file of the data we want to analyse
-- data is the phenotype file, likely a superset of PLINK data
-- cov_list : list of comma-separated strings
-- pheno : single phenotype label
-- dataout : name of a phenotype file that gets created by the script. There is guaranteed no missing phenotype data in the file, and
   only contains individuals described in the fam file -- not guaranteed that it's in the same order as the fam file
-- ind_list: subset of the individuals in the fam file for (a) no missing genotype data (b) no missing phenotype data

"""

import sys
import pandas as pd
import argparse
import numpy as np




def check_missing(x):
    if x in ["NA","na","null","-9"]  or ((type(x)==type("x")) and (len(x)<1)) or x==-9 or x==r"\N" or x=="-":
        return 1
    else:
        return 0

def parseArguments():
    parser = argparse.ArgumentParser(description='fill in missing bim values')
    parser.add_argument('--inp_fam',type=str,required=True)
    parser.add_argument('--data',type=str,required=True,help="File with phenotype and covariate data")
    parser.add_argument('--cov_list', type=str,help="comma separated list of covariates",default="")
    parser.add_argument('--pheno',type=str,required=True,help="single phenotype")
    parser.add_argument('--dataout',type=str,required=True,help="File with output pheno")
    parser.add_argument('--lindout',type=str,required=True,help="File with list of good individuals FID IID ")
    args = parser.parse_args()
    return args

def getColNames(label_str):
   col_names = []
   col_fns   = []
   for lab in label_str:
       det = lab.split("/")
       if len(det)>1:
           col_names.append(det[0])
           col_fns.append(eval(det[1]))
       else:
           col_names.append(lab)
           col_fns.append(False)
   return col_names,col_fns

def errorMessage10(phe):
    print("""

    A problem has been detected in file <%s> column <%s>.

    There is some invalid data. I regret I can't tell you which row.


    Please check -- the data should be numeric only.


    If there is missing data, please use   NA



    """%(args.data,phe))


args = parseArguments()

TAB =chr(9)
EOL=chr(10)


phenos     = [args.pheno]
use        = ["FID","IID"]
if args.cov_list:
    phenos+= list(set([x for x in args.cov_list.split(",") if len(x)>0]))

#read of fam 
readfam=open(args.inp_fam)
listeFID=[]
for Lines in readfam :
    SplL=Lines.split()
    listeFID.append(SplL[0]+"-"+SplL[1])
readfam.close()

readdata=open(args.data)
EnteteFile=readdata.readline().split()
UsePos=[]
## search position of interest
useAll=use+phenos
for Ent in useAll :
    if Ent not in EnteteFile :
      sys.exit(EOL*3+"column %s not find in file %s"%(Ent,args.data))
    UsePos.append(EnteteFile.index(Ent))

writenewdata=open(args.dataout,'w')
writenewdata.write("\t".join(useAll)+"\n")
writelind=open(args.lindout,'w')
for line in readdata : 
    splline=line.split()    
    if splline[UsePos[0]]+"-"+splline[UsePos[1]] in listeFID:
       if sum([check_missing(splline[x]) for x in UsePos[2:]])==0 :
          NewInf=[splline[x] for x in UsePos]
          writelind.write(NewInf[0]+"\t"+NewInf[1]+"\n")
          writenewdata.write("\t".join(NewInf)+"\n")

writenewdata.close()
writelind.close()
readdata.close()







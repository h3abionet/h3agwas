#!/usr/bin/env python3

# Converts PLINK covariate and fam file into a covariate file for Gemma

import sys
import pandas as pd
import argparse
import numpy as np

def parseArguments():
    parser = argparse.ArgumentParser(description='fill in missing bim values')
    parser.add_argument('--inp_fam',type=str,required=True)
    parser.add_argument('--data',type=str,required=True,help="File with phenotype and covariate data")
    parser.add_argument('--cov_list', type=str,required=True,help="comma separated list of covariates")
    parser.add_argument('--pheno',type=str,required=True,help="comma separated list of  pheno column")
    parser.add_argument('--phe_out', type=str,help="output fam file")
    parser.add_argument('--cov_out', type=str,help="output covariate file")
    args = parser.parse_args()
    return args


def check_missing(x):
    if x in ["NA","na","null","-9"]  or ((type(x)==type("x")) and (len(x)<1)) or x==-9 or x==r"\N" or x=="-":
        return "NA"
    else:
        return x

args = parseArguments()

TAB =chr(9)

covariates = args.cov_list.split(",")
use = ["FID","IID"]+covariates

phenos_0 = args.pheno.split(",")
phenos = []
for lab in phenos_0:
    det = lab.split("/")
    if len(det)>1:
        phenos.append(det[0])
    else:
        phenos.append(lab)


famd  = pd.read_csv(args.inp_fam,header=None,delim_whitespace=True,names=["FID","IID","FAT","MAT","SEX","CC"])


datad = pd.read_csv(args.data,delim_whitespace=True,usecols=phenos+use)
columns=datad.columns()
if "FID" not in columns or "IID" not in columns:
    sys.exit(EOL*3+"It is mandatory for FID IID to be columns of the file %s"%args.data+EOL*3)


for phe in phenos_0:
    det = lab.split("/")
    phe = det[0]
    if phe not in columns:
        sys.exit((EOL*3+"<%s> given as phenotype, but is not a column of the file <%s>"+EOL+EOL)%(phe,args.datad))
    if len(det)>1:
        fn = eval(det[1])
        datad[phe]=fn(datad[phe])
    datad[phe]=datad[phe].apply(check_missing)

merge = pd.merge(famd,datad,how="left",on=["FID","IID"])
merge["intercept"]=1



merge.reindex(["FID","IID","intercept"]+covariates)
merge.to_csv(args.cov_out,sep=TAB,columns=["intercept"]+covariates,header=False,index=False,na_rep=-9)

merge.to_csv(args.phe_out,sep=TAB,columns=phenos,header=False,index=False,na_rep=-9)

print(" ".join(list(map (lambda a: str(a[0]+1)+"-"+a[1], enumerate(phenos) ))),end="")

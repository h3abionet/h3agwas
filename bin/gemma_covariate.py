#!/usr/bin/env python3

# Converts PLINK covariate and fam file into a covariate file for Gemma

import sys
import pandas as pd
import argparse
import numpy as np


EOL=chr(10)

def parseArguments():
    parser = argparse.ArgumentParser(description='fill in missing bim values')
    parser.add_argument('--inp_fam',type=str,required=True)
    parser.add_argument('--data',type=str,required=True,help="File with phenotype and covariate data")
    parser.add_argument('--cov_list', type=str,help="comma separated list of covariates",default="")
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

phenos     = args.pheno.split(",")
use        = ["FID","IID"]
if args.cov_list:
    covariates = args.cov_list.split(",")
    use = use+covariates
else:
    covariates = []
    


pheno_labels, pheno_transform  = getColNames(phenos)
covar_labels, cover_transform  = getColNames(use)

usecols = covar_labels+pheno_labels

datad = pd.read_csv(args.data,delim_whitespace=True,usecols=usecols)
columns = datad.columns

if "FID" not in columns or "IID" not in columns:
    sys.exit(EOL*3+"It is mandatory for FID IID to be columns of the file %s"%args.data+EOL*3)

for (label, transform) in zip(pheno_labels+covar_labels, pheno_transform+cover_transform) :
    if label not in columns:
        sys.exit((EOL*3+"<%s> given as label, but is not a column of the file <%s>"+EOL+EOL)%(label,args.datad))
    if transform:
        try:
            datad[label]=transform(datad[label])
        except:
            errorMessage10(label)
            sys.exit(10)
    datad[label]=datad[label].apply(check_missing)


famd  = pd.read_csv(args.inp_fam,header=None,delim_whitespace=True,names=["FID","IID","FAT","MAT","SEXFAM","CC"])


merge = pd.merge(famd,datad,how="left",on=["FID","IID"])
merge["intercept"]=1
merge.reindex(["FID","IID","intercept"]+covariates)
merge.to_csv(args.cov_out,sep=TAB,columns=["intercept"]+covariates,header=False,index=False,na_rep="NA")
merge.to_csv(args.phe_out,sep=TAB,columns=pheno_labels,header=False,index=False,na_rep="NA")
print(" ".join(list(map (lambda a: str(a[0]+1)+"-"+a[1], enumerate(phenos) ))),end="")

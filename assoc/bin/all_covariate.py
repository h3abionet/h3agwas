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
    parser.add_argument('--lind',type=str,required=False)
    parser.add_argument('--data',type=str,required=True,help="File with phenotype and covariate data")
    parser.add_argument('--cov_list', type=str,help="comma separated list of covariates",default="")
    parser.add_argument('--cov_type', type=str,help="comma separated list of covariates",default="")
    parser.add_argument('--pheno',type=str,required=True,help="comma separated list of  pheno column")
    parser.add_argument('--phe_out', type=str,help="output fam file")
    parser.add_argument('--cov_out', type=str,help="output covariate file")
    parser.add_argument('--covqual_file', type=str,help="output covariate file")
    parser.add_argument('--cov_file', type=str,help="output covariate file")
    parser.add_argument('--gxe_out', type=str,help="output gxe file (gemma use)")
    parser.add_argument('--gxe', type=str,help="gxe covariate (gemma use)")
    parser.add_argument('--form_out', type=int,help="format output : 1:Gemma, 2:boltlmm, 3:FastLmm, 4:gcta", required=True)
    args = parser.parse_args()
    return args


def check_missing(x, MissingOut):
    if x in ["NA","na","null","-9"]  or ((type(x)==type("x")) and (len(x)<1)) or x==-9 or x==r"\N" or x=="-":
        return MissingOut
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
use= ["FID","IID"]
useI        = ["FID","IID"]
if args.cov_list:
    covariates = args.cov_list.split(",")
    use = useI+covariates
else:
    covariates = []

if args.gxe:
   gxe=[args.gxe]
   gxe_use = useI+gxe
else :
   gxe=[]
  

if args.form_out==1 :
   MissingOut="NA"
elif  args.form_out==2:
    MissingOut="NA" 
elif  args.form_out==3:
    MissingOut="-9" 
elif  args.form_out==4:
    MissingOut="NA" 
else :
    print("--form_out : "+str(args.form_out)+" not define")
    sys.exit(11)
    


pheno_labels, pheno_transform  = getColNames(phenos)
covar_labels, cover_transform  = getColNames(use)
gxe_labels, gxe_transform  = getColNames(gxe)

usecols = covar_labels+pheno_labels+gxe_labels

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
    datad[label]=datad[label].apply(check_missing, MissingOut=MissingOut)

famd  = pd.read_csv(args.inp_fam,header=None,delim_whitespace=True,names=["FID","IID","FAT","MAT","SEXFAM","CC"])
#if args.lind :
#   indtokepp=pd.read_csv(args.lind,delim_whitespace=True,header=None, names=["FID","IID"])
#   famd=pd.merge(famd, indtokepp,how="inner",on=["FID","IID"])
merge = pd.merge(famd,datad,how="left",suffixes=["_f",""],on=["FID","IID"])
# for gemma
if args.form_out == 1 : 
   merge["intercept"]=1
   merge.reindex(["FID","IID","intercept"]+covariates)
   merge.to_csv(args.cov_out,sep=TAB,columns=["intercept"]+covariates,header=False,index=False,na_rep=MissingOut)
   merge.to_csv(args.phe_out,sep=TAB,columns=pheno_labels,header=False,index=False,na_rep=MissingOut)
   merge.reindex(["FID","IID"]+gxe)
   merge.to_csv(args.gxe_out,sep=TAB,columns=gxe_labels,header=False,index=False,na_rep=MissingOut)
elif  args.form_out == 2 :
   merge.reindex(["FID","IID"])
   merge.to_csv(args.phe_out,sep=TAB,columns=["FID","IID"]+pheno_labels+covariates,header=True,index=False,na_rep=MissingOut)
elif  args.form_out == 3 :
   merge.reindex(["FID","IID"])
   merge.to_csv(args.phe_out,sep=TAB,columns=["FID","IID"]+pheno_labels,header=False,index=False,na_rep=MissingOut)
   merge.to_csv(args.cov_out,sep=TAB,columns=["FID","IID"]+covariates,header=False,index=False,na_rep=MissingOut)
elif args.form_out == 4 :
   merge.reindex(["FID","IID"])
   merge.to_csv(args.phe_out,sep=TAB,columns=["FID","IID"]+pheno_labels,header=True,index=False,na_rep=MissingOut)
   if not args.cov_type :
      merge.to_csv(args.cov_out,sep=TAB,columns=["FID","IID"]+covariates,header=True,index=False,na_rep=MissingOut)
   else :
      splcov=args.cov_type.split(',')      
      covtmp=covariates
      if len(splcov)!=len(splcov):
         print('error cov '+ covariates+' have different len than type cov'+ args.cov_type)
         sys.exit(10)
      covquant=[]
      covqual=[]
      cmt=0
      for typecov in splcov:
         if typecov =="1":
             covquant.append(covtmp[cmt])
         elif typecov== "0":
             covqual.append(covtmp[cmt])
         else :
             print('error : type covariable different of 0 (qualitatif) and 1 (quantitatif) '+typecov)
         cmt+=1
      if len(covquant)>0:
         merge.to_csv(args.cov_out,sep=TAB,columns=["FID","IID"]+covquant,header=True,index=False,na_rep=MissingOut)
      merge.to_csv(args.covqual_file,sep=TAB,columns=["FID","IID"]+covqual,header=True,index=False,na_rep=MissingOut)


print(" ".join(list(map (lambda a: str(a[0]+1)+"@@@"+a[1], enumerate(phenos) ))),end="")

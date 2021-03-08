#!/usr/bin/env python3
import pandas as pd
import argparse
import numpy as np
import sys

def is_missing(x):
    if x in ["NA","na","null","-9"]  or ((type(x)==type("x")) and (len(x)<1)) or x==-9 or x==r"\N" or x=="-":
        return True
    else:
        return False


''' 
normalise phenotype by some variable
'''

TAB =chr(9)
def normalise(x, minx, maxx, newmin, newmax, na_out):
    if is_missing(x) :
      return na_out
    return (float(x)-minx)/(maxx-minx) * (newmax-newmin) + newmin


def parseArguments():
    parser = argparse.ArgumentParser(description='Normalise variable in function of data and parameters, if not --cov_info, merge file together data and pheno')
    parser.add_argument('--data',type=str, required=True, help="phenotype files contains pheno initial, and covariable")
    parser.add_argument("--data_sim", required=True, help="phenotype simulate with FID and IID")
    parser.add_argument('--cov_info',type=str, required=False, help="info of covariable cov1=x,cov2=y...")
    parser.add_argument('--out',type=str, required=True, help="file out data")
    parser.add_argument('--rangenorm',type=str, required=False, help="values to normarlise (use min,max of residual)")
    parser.add_argument('--intercept',type=float, required=False, help="intercept used for normalisation", default=0)
    parser.add_argument('--na_out',type=str, required=True, help="na value used to output")
    args = parser.parse_args()
    return args

args=parseArguments()
data  = pd.read_csv(args.data,delim_whitespace=True)
data_phenosim=pd.read_csv(args.data_sim, delim_whitespace=True)
na_out=args.na_out
if args.cov_info==None :
   merge_pheno=pd.merge(data,data_phenosim,how="inner",on=["FID","IID"])
else :
   RangeNorm=[float(x) for x in args.rangenorm.split(",")]
   listcov=[x.split("=")[0]   for x in args.cov_info.split(",")]
   valcov=[float(x.split("=")[1])   for x in args.cov_info.split(",")]
   ColNames=["FID","IID"]
   ColNames+=valcov
   KeysPheno=data_phenosim.keys()[2]
   minpheno=min(data_phenosim[KeysPheno])
   maxpheno=max(data_phenosim[KeysPheno])
   NewPhen=KeysPheno+"_N"
   data_phenosim[NewPhen]=data_phenosim[KeysPheno].apply(normalise, minx=minpheno,maxx=maxpheno, newmin=RangeNorm[0],newmax=RangeNorm[1], na_out=na_out)
   merge_pheno=pd.merge(data,data_phenosim,how="inner",on=["FID","IID"])
   Inter=args.intercept
   ## merge by list
   finalval=[]
   for index,row in merge_pheno.iterrows():
      valPhenoI=row[NewPhen]
      if is_missing(valPhenoI) :
         valPhenoI=na_out
      else :
        for cmt in range(len(listcov)) :
           if is_missing(valcov[cmt]) or valPhenoI==na_out: 
              valPhenoI=na_out
           else :
              valPhenoI+=row[listcov[cmt]]*valcov[cmt]
        if is_missing(valcov[cmt])==False:
           valPhenoI+=Inter
      finalval.append(valPhenoI)
   merge_pheno[KeysPheno]=finalval
merge_pheno.to_csv(args.out,sep=TAB,header=True,index=False, na_rep=na_out)




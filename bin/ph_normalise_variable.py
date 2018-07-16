#!/usr/bin/env python3
import pandas as pd
import argparse
import numpy as np
import sys


''' 
normalise phenotype by some variable
'''

TAB =chr(9)
def normalise(x, minx, maxx, newmin, newmax):
    return (float(x)-minx)/(maxx-minx) * (newmax-newmin) + newmin


def parseArguments():
    parser = argparse.ArgumentParser(description='Normalise variable in function of data and parameters, if not --cov_info, merge file together data and pheno')
    parser.add_argument('--data',type=str, required=True, help="phenotype files contains pheno initial, and covariable")
    parser.add_argument("--data_sim", required=True, help="")
    parser.add_argument('--cov_info',type=str, required=False, help="info of covariable cov1=x,cov2=y...")
    parser.add_argument('--out',type=str, required=True, help="file out data")
    parser.add_argument('--rangenorm',type=str, required=False, help="values to normarlise (use min,max of residual)")
    parser.add_argument('--intercept',type=float, required=False, help="values to normarlise (use min,max of residual)", default=0)
    args = parser.parse_args()
    return args

args=parseArguments()
data  = pd.read_csv(args.data,delim_whitespace=True)
data_phenosim=pd.read_csv(args.data_sim, delim_whitespace=True)
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
   data_phenosim[NewPhen]=data_phenosim[KeysPheno].apply(normalise, minx=minpheno,maxx=maxpheno, newmin=RangeNorm[0],newmax=RangeNorm[1])
   merge_pheno=pd.merge(data,data_phenosim,how="inner",on=["FID","IID"])
   Inter=args.intercept
   ## merge by list
   finalval=[]
   for index,row in merge_pheno.iterrows():
      valPhenoI=row[NewPhen]
      for cmt in range(len(listcov)) :
         valPhenoI+=row[listcov[cmt]]*valcov[cmt]
      finalval.append(valPhenoI+Inter)
   merge_pheno[KeysPheno]=finalval
merge_pheno.to_csv(args.out,sep=TAB,header=True,index=False)


        






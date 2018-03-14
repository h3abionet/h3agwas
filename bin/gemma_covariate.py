#!/usr/bin/env python3

# Converts PLINK covariate and fam file into a covariate file for Gemma

import sys
import pandas as pd
import argparse

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
    if len(x)<1 or x=="-9" or x==r"\N" or x=="-":
        return "NA"
    else:
        return x

args = parseArguments()

TAB =chr(9)

covariates = args.cov_list.split(",")
use = ["FID","IID"]+covariates

phenos = args.pheno.split(",")

famd  = pd.read_csv(args.inp_fam,header=None,delim_whitespace=True,names=["FID","IID","FAT","MAT","SEX","CC"])

datad = pd.read_csv(args.data,delim_whitespace=True,usecols=phenos+use)
for phe in phenos:
    datad[phe]=datad[phe].apply(check_missing)

merge = pd.merge(famd,datad,how="left",on=["FID","IID"])
merge["intercept"]=1



merge.reindex(["FID","IID","intercept"]+covariates)
merge.to_csv(args.cov_out,sep=TAB,columns=["intercept"]+covariates,header=False,index=False)

merge.to_csv(args.phe_out,sep=TAB,columns=phenos,header=False,index=False)


phe_cols = [col for col in merge.columns if col in phenos]
print(" ".join(list(map (lambda a: str(a[0]), enumerate(phe_cols) ))))

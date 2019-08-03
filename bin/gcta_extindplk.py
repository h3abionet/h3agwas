#!/usr/bin/env python3
''' 
extract list of individu in function of phenotype with data file
'''

from  subprocess import CalledProcessError
import subprocess
import os
import argparse
import numpy as np
import pandas as pd
import sys

def parseArguments():
    parser = argparse.ArgumentParser(description='format file for gcta, append N and frequencie if not present using bed file')
    parser.add_argument('--data',type=str,required=True, help="phenotypes files")
    parser.add_argument('--out', type=str,help="out of tex file")
    parser.add_argument('--pheno',type=str,required=False,help="phenotype header")
    args = parser.parse_args()
    return args

args = parseArguments()

data=pd.read_csv(args.data,delim_whitespace=True)
if args.pheno :
   data=data[pd.isna(data[args.pheno])==False]
data=data.iloc[:,0:2] 
data.to_csv(args.out,sep="\t",header=False,index=False,na_rep="NA")





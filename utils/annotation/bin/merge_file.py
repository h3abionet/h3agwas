#!/usr/bin/env python3

import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import os

TAB =chr(9)

def parseArguments():
    parser = argparse.ArgumentParser(description='Produces Manhatten, QQ plot and supporting tex file  for output from some tools')
    parser.add_argument('--list_pos',type=str,required=True, help="contains pos to extract (2 column)")
    parser.add_argument('--list_file',type=str,required=True, help="annovar format file")
    parser.add_argument('--out', type=str,help="out of tex file",required=True)
    args = parser.parse_args()
    return args
args = parseArguments()



DataI=pd.read_csv(args.list_pos,delim_whitespace=True,header=None)
DataI.columns = ["#Chr","Start","End","Ref","Alt"]
DataI.drop(["End"])
LireFileAnno=args.list_file.split(',') #open(args.list_file).readlines()

for File in LireFileAnno :
    dataAnnot=pd.read_csv(File.replace("\n",""),delim_whitespace=True)
    if 'End' in dataAnnot.columns :
      del dataAnnot['End']
    DataI = pd.merge(DataI,dataAnnot,how="left",on=["#Chr","Start","Ref","Alt"])

DataI.to_csv(args.out,sep=TAB,header=True,index=False, na_rep="NA")

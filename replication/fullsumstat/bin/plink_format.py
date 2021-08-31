#!/usr/bin/env python3
''' 
format file for gcta, append freq and n if need and bfile
'''
from  subprocess import CalledProcessError
import subprocess
import os
import argparse
import numpy as np
import pandas as pd
import sys

EOL = chr(10)
#SNP A1 A2 freq b se p N 
def parseArguments():
    parser = argparse.ArgumentParser(description='fill in missing bim values')
    parser.add_argument('--inp_asso',type=str,required=True, help="association files")
    parser.add_argument('--out', type=str,help="out of tex file",default="test.tex")
    parser.add_argument('--rs_header',type=str,required=True,help="rs header in inp files")
    parser.add_argument('--pval_header',type=str,required=True,help="pvalue header in inp files")
    parser.add_argument('--freq_header',type=str,required=False,help="freq header in inp files",default=None)
    parser.add_argument('--a1_header',type=str,required=True,help="a1 header in inp files")
    parser.add_argument('--a2_header',type=str,required=True,help="a2 header in inp files")
    parser.add_argument('--se_header',type=str,required=True,help="se header in inp files")
    parser.add_argument('--n_header',type=str,help="n header in inp files", default=None)
    parser.add_argument('--chr',type=str,help="n header in inp files", default=None)
    parser.add_argument('--chro_header',type=str,required=True,help="n header in inp files")
    parser.add_argument('--bp_header',type=str,help="bp header")
    parser.add_argument('--beta_header',type=str,required=True,help="beta header in inp files")
    parser.add_argument('--chr_pos',type=str,required=False,help="beta header in inp files")
    args = parser.parse_args()
    return args

args = parseArguments()
inp     = args.inp_asso

rs_head=args.rs_header
beta_head=args.beta_header
freq_head=args.freq_header
pval_head=args.pval_header
se_head=args.se_header
a1_head=args.a1_header
a2_head=args.a2_header
se_head=args.se_header
n_head=args.n_header
IsFreq=False
IsN=False
if freq_head :
  IsFreq=True
if n_head :
  IsN=True

## first step : if
result = pd.read_csv(inp,delim_whitespace=True)
result[args.chro_header] = result[args.chro_header].astype(str)

if args.chr :
   result=result[result[args.chro_header]==args.chr]

if args.chr_pos :
   chrpos= pd.read_csv(args.chr_pos,delim_whitespace=True,header=None)
   chrpos.columns=[args.chro_header, args.bp_header]
   chrpos[args.chro_header] = chrpos[args.chro_header].astype(str)
   result=pd.merge(chrpos,  result, on=[args.chro_header, args.bp_header])
   


#Head=[rs_head, a1_head,a2_head,beta_head, se_head,pval_head]
#NewHead=["SNP","A1","A2","b","se","p"]
NewHead=["SNP"]
Head=[rs_head]
NewHead.append("CHR")
Head.append(args.chro_header)

NewHead.append("BP")
Head.append(args.bp_header)

NewHead+=["A1","A2"]
Head+=[a1_head,a2_head]

NewHead+=["BETA","SE"]
Head+=[beta_head, se_head]
if pval_head :
   NewHead+=["P"]
   Head+=[pval_head]

if IsN :
   Head.append(n_head)
   NewHead.append("N")

result=result.loc[:,Head]
result.columns=NewHead


result.to_csv(args.out,sep=" ",header=True,index=False,na_rep="NA")



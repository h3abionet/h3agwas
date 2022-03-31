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
    parser.add_argument('--typegwas',type=str,required=False,help="rs header in inp files")
    parser.add_argument('--rs_header',type=str,required=False,help="rs header in inp files")
    parser.add_argument('--pval_header',type=str,required=False,help="pvalue header in inp files")
    parser.add_argument('--freq_header',type=str,required=False,help="freq header in inp files",default=None)
    parser.add_argument('--a1_header',type=str,required=False,help="a1 header in inp files")
    parser.add_argument('--a2_header',type=str,required=False,help="a2 header in inp files")
    parser.add_argument('--se_header',type=str,required=False,help="se header in inp files")
    parser.add_argument('--n_header',type=str,help="n header in inp files", default=None)
    parser.add_argument('--chr',type=str,help="n header in inp files", default=None)
    parser.add_argument('--chro_header',type=str,required=False,help="n header in inp files")
    parser.add_argument('--bp_header',type=str,help="bp header")
    parser.add_argument('--beta_header',type=str,required=False,help="beta header in inp files")
    parser.add_argument('--chr_pos',type=str,required=False,help="beta header in inp files")
    parser.add_argument('--bim',type=str,required=False,help="beta header in inp files")
    args = parser.parse_args()
    return args

args = parseArguments()
inp     = args.inp_asso

#headchr=chr;headpval="p_wald";headfreq=af;headA1=allele1;headA2=allele2;headrs=rs;headbeta="beta";headbp="ps";headse="se"
if args.typegwas=='gemma'or args.typegwas=="gemmaN":
  rs_head='rs'
  chro_head='chr'
  bp_head='ps'
  beta_head='beta'
  freq_head='af'
  pval_head='p_wald'
  se_head='se'
  a1_head='allele1'
  a2_head='allele0'
  n_head='N'
elif args.typegwas=='metal' :
  rs_head='SNP'
  chro_head='CHRO'
  bp_head='BP'
  beta_head=None
  freq_head='Freq1'
  pval_head='P-value'
  se_head=None
  a1_head='Allele1'
  a2_head='Allele2'
  n_head='Weight'
  ##SNP    CHRO    BP      MarkerName      Allele1 Allele2 Freq1   FreqSE  MinFreq MaxFreq Weight  Zscore  N       P-value Direction
else :
  chro_head=args.chro_header
  rs_head=args.rs_header
  beta_head=args.beta_header
  freq_head=args.freq_header
  pval_head=args.pval_header
  se_head=args.se_header
  a1_head=args.a1_header
  a2_head=args.a2_header
  se_head=args.se_header
  n_head=args.n_header
  bp_head=args.bp_header


IsFreq=False
IsN=False
if freq_head :
  IsFreq=True
if n_head :
  IsN=True

## first step : if
result = pd.read_csv(inp,delim_whitespace=True)
result[chro_head] = result[chro_head].astype(str)

if args.chr :
   result=result[result[chro_head]==args.chr]

if args.chr_pos :
   chrpos= pd.read_csv(chr_pos,delim_whitespace=True,header=None)
   chrpos.columns=[chro_head, bp_head]
   chrpos[chro_head] = chrpos[chro_head].astype(str)
   result=pd.merge(chrpos,  result, on=[chro_head, bp_head])
   

result[a1_head] = result[a1_head].str.upper()
result[a2_head] = result[a2_head].str.upper()

#Head=[rs_head, a1_head,a2_head,beta_head, se_head,pval_head]
#NewHead=["SNP","A1","A2","b","se","p"]
NewHead=["SNP2"]
Head=[rs_head]
NewHead.append("CHR")
Head.append(chro_head)

NewHead.append("BP")
Head.append(bp_head)

NewHead+=["A1","A2"]
Head+=[a1_head,a2_head]

if beta_head :
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

bim=pd.read_csv(args.bim,delim_whitespace=True,header=None)
bim.columns= ["CHR", "SNP", "CM", "BP", "A1","A2"]
bim=bim[['CHR', 'SNP', 'BP']]

bim['CHR']=bim['CHR'].astype(str)
result['CHR']=result['CHR'].astype(str)

result=bim.merge(result, left_on=['CHR','BP'], right_on=['CHR','BP'])



result.to_csv(args.out,sep=" ",header=True,index=False,na_rep="NA")



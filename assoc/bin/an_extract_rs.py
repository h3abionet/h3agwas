#!/usr/bin/env python3

# Converts PLINK covariate and fam file into a covariate file for Gemma

import sys
import pandas as pd
import argparse
import numpy as np

EOL=chr(10)

def parseArguments():
    parser = argparse.ArgumentParser(description='extract rs in gwas file')
    parser.add_argument('--inp_resgwas',type=str,required=True)
    parser.add_argument('--inp_rs',type=str,required=False, help="file contains rs", default="")
    parser.add_argument('--list_rs',type=str,required=False, help="list rs", default="")
    parser.add_argument('--chro_header',type=str,required=True,help="chro header in inp files")
    parser.add_argument('--pos_header',type=str,required=True,help="pos header in inp files")
    parser.add_argument('--rs_header',type=str,required=True,help="rs header in inp files")
    parser.add_argument('--pval_header',type=str,required=True,help="pvalue header in inp files")
    parser.add_argument('--beta_header',type=str,required=False,help="beta header in inp files")
    parser.add_argument('--freq_header',type=str,required=False,help="frequencies header in inp files")
    parser.add_argument('--around_rs',type=float,required=True,help="around rs (pb)")
    parser.add_argument('--maf',type=float,default=0.0,help="minor allele frequencies")
    parser.add_argument('--out_head',type=str,default="out",help="around rs (pb)")
    args = parser.parse_args()
    return args

def read_rs(FileRs) :
    ListeRs=[]
    read=open(FileRs)
    for x in read :
        ListeRs.append(x.replace('\n', '').split()[0])
    return ListeRs

args = parseArguments()

if args.inp_rs != "" :
   list_rs=read_rs(args.inp_rs)
elif args.list_rs != "" :
     list_rs=args.list_rs.split(",")
else :
      sys.exit("args : inp_rs and list_rs not found\nexit\n")

result = pd.read_csv(args.inp_resgwas,delim_whitespace=True, dtype={args.chro_header:str})
#result=result[(result[args.freq_header]>args.maf) & (result[args.freq_header]<1-args.maf)]

sub_result=result.loc[result[args.rs_header].isin(list_rs)]
TAB=chr(9)
if args.freq_header :
   PosCol=[args.chro_header, args.pos_header,  args.pos_header, args.rs_header, args.pval_header, args.freq_header]
else :
   PosCol=[args.chro_header, args.pos_header,  args.pos_header, args.rs_header, args.pval_header]
for x in sub_result[args.rs_header] :
    out_file=args.out_head+"_around.stat"
    out_gwas=args.out_head+"_gwas.stat"
    out_info=args.out_head+"_info.stat"
    infors=sub_result[sub_result[args.rs_header]==x]
    infors.to_csv(out_gwas, sep=TAB, header=True, index=False,na_rep="NA")
    chro=infors[args.chro_header].tolist()[0]
    pos=infors[args.pos_header].tolist()[0]
    posbegin=pos-args.around_rs
    posend=pos+args.around_rs
    small=result[(result[args.chro_header]==chro) & (result[args.pos_header]>=posbegin) & (result[args.pos_header]<=posend)]
    maf=args.maf
    if args.freq_header :
       small=small[(small[args.rs_header]==x) | ((small[args.freq_header]>maf) & (small[args.freq_header]<(1-maf)))]
    # chrom, start, end, marker ID, and p-value 
    verysmall=infors[[args.rs_header,args.chro_header, args.pos_header]]
    verysmall.to_csv(out_info, sep=TAB, header=False, index=False,na_rep="NA")
    small=small[PosCol] 
    if args.freq_header :
       small.columns=["#CHROM","BEGIN","END","MARKER_ID","PVALUE","MAF"]
    else :
       small.columns=["#CHROM","BEGIN","END","MARKER_ID","PVALUE"]
    small.BEGIN = small.BEGIN.astype(int)
    small.END = small.END.astype(int)
    #small['#CHROM'] = small['#CHROM'].astype(int)
    small.to_csv(out_file, sep=TAB, header=True, index=False,na_rep="NA")


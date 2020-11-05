#!/usr/bin/env python3

import sys
import os
import argparse
import math
import numpy as np
import pandas as pd
def addedplkinfo(args, datai, ChroHead, BpHead):
   bim=pd.read_csv(args.bfile+".bim", header=None,delim_whitespace=True)
   #1	1:725499	0	725499	T	A
   bim.columns=['Chro','rsold','O0','Bphead','O1','02']
   bim['Chro']=bim['Chro'].astype('str')
   if args.chro :
      bim=bim[bim['Chro']==args.chro]
   bim=bim.iloc[:,[0,1,3]]
   bim.columns=[ChroHead,'rsold', BpHead]
   plkfreqfil=os.path.basename(args.bfile)
   if args.bfile==None :
     print("no header for n or freq and bfile")
     sys.exit()
   Cmd=args.bin_plk+" -bfile "+args.bfile+" --freq --keep-allele-order "
   if args.chro :
     Cmd+=" --chr "+args.chro
     plkfreqfil=plkfreqfil+"_"+args.chro
   if args.keep :
     Cmd+=" --keep "+args.keep
   Cmd+=" --threads "+str(args.threads)
   Cmd+=" --out "+plkfreqfil
   os.system(Cmd)
   data_n=pd.read_csv(plkfreqfil+".frq",delim_whitespace=True)
   del data_n['CHR']
   data_n=data_n.merge(bim, left_on='SNP', right_on='rsold')
   data_n['N']=data_n['NCHROBS']/2
   if args.N_head==None and args.freq_head==None:
      data_n=data_n[[ChroHead, BpHead,"N",'MAF']]
      data_n.columns=[ChroHead, BpHead, args.Nnew_head, args.freqnew_head]
   elif args.N_head==None :
      data_n=data_n[[ChroHead, BpHead,"N"]]
      data_n.columns=[ChroHead, BpHead, args.Nnew_head]
   elif args.freq_head==None:
      data_n=data_n[[ChroHead, BpHead,'MAF']]
      data_n.columns=[ChroHead, BpHead, args.freqnew_head]
   datai=datai.merge(data_n, on=[ChroHead, BpHead],how='left')
   return datai




def parseArguments():
    parser = argparse.ArgumentParser(description='transform file and header')
    parser.add_argument('--input_gwas',type=str,required=True, help="input gwas")
    parser.add_argument('--input_rs',type=str,required=True, help="input new rs")
    parser.add_argument('--out_file',type=str,required=True,help="output file")
    parser.add_argument('--chro',type=str,required=False,help="specific chromosome")
    parser.add_argument('--chro_head',type=str,required=True,help="header of chromosome")
    parser.add_argument('--bp_head',type=str,required=True,help="head of postion")
    parser.add_argument('--rs_head',type=str,required=True,help="new head of rs")
    parser.add_argument('--a1_head',type=str,required=True,help="head of postion")
    parser.add_argument('--a2_head',type=str,required=True,help="new head of rs")
    parser.add_argument('--freq_head',type=str,required=False,help="new head of rs")
    parser.add_argument('--freqnew_head',type=str,default='Freq',required=False,help="new head of rs")
    parser.add_argument('--N_head',type=str,required=False,help="")
    parser.add_argument('--Nnew_head',type=str,default='N',required=False,help="")
    parser.add_argument('--bfile',type=str,required=False,help="plink file")
    parser.add_argument('--keep',type=str,required=False,help="option keep for plink")
    parser.add_argument('--threads',type=str,default=1,required=False,help="option threads for plink")
    parser.add_argument('--N_value', type=str,required=False,help="added a N values")
    parser.add_argument('--bin_plk',type=str,required=False,help="plink binary", default='plink')
    args = parser.parse_args()
    return args



args=parseArguments()
ChroHead=args.chro_head
BPHead=args.bp_head
RsHead=args.rs_head
gwasres = pd.read_csv(args.input_gwas,delim_whitespace=True)
inforrs=pd.read_csv(args.input_rs,header=None, delim_whitespace=True)
inforrs.columns =[ChroHead,BPHead,args.a1_head+'_tmp', args.a2_head+'_tmp', RsHead, "a1_old"]
gwasres=inforrs.merge(gwasres, how='right',on=[ChroHead,BPHead])
balise=gwasres[RsHead].isnull() | ((((gwasres[args.a1_head]==gwasres[args.a1_head+'_tmp']) & (gwasres[args.a2_head]==gwasres[args.a2_head+'_tmp'])) | ((gwasres[args.a2_head]==gwasres[args.a1_head+'_tmp']) & (gwasres[args.a1_head]==gwasres[args.a2_head+'_tmp'])))==False)
gwasres.loc[balise,RsHead]=gwasres.loc[balise,ChroHead].astype('str')+":"+gwasres.loc[balise,BPHead].astype('str')
gwasres[ChroHead]=gwasres[ChroHead].astype('str')
if args.bfile and (args.N_head==None or args.freq_head==None):
   gwasres=addedplkinfo(args, gwasres, ChroHead,BPHead) 
elif args.N_value and args.N_head==None:
  gwasres[args.Nnew_head]=args.N_value
gwasres.to_csv(args.out_file,sep=" ",header=True,index=False,na_rep="NA")



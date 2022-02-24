#!/usr/bin/env python3
import random
import os
import math

''' 
Produces Manhatten, QQ plot and supporting tex file  for output from some tools
'''

import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

def parseArguments():
    parser = argparse.ArgumentParser(description='Produces Manhatten, QQ plot and supporting tex file  for output from some tools')
    parser.add_argument('--file_gxe',type=str,required=True, help="association files")
    parser.add_argument('--bfile',type=str,required=True,help="pheno names")
    parser.add_argument('--pheno_file',type=str,help="pheno names")
    parser.add_argument('--pheno',type=str,help="pheno names")
    parser.add_argument('--pheno_gxe',type=str,help="pheno names")
    parser.add_argument('--freq1',type=str,help="pheno names")
    parser.add_argument('--freq2',type=str,help="pheno names")
    parser.add_argument('--freqAll',type=str,help="pheno names")
    parser.add_argument('--out', type=str,help="",default="out")
    parser.add_argument('--plk_cores', type=str,help="",default="2")
    parser.add_argument('--gwas_chr', type=str,help="",default='chr')
    parser.add_argument('--gwas_ps', type=str,help="",default='ps')
    parser.add_argument('--gwas_rs', type=str,help="",default='rs')
    parser.add_argument('--dir_tmp', type=str, help="")
    args = parser.parse_args()
    return args

Test=False

#args={'file_gxe':'All/all_imputed-meancIMT_Res.qassoc.final.gxe.modif','FilefreqAll':'All/Freq/CIMT_All.pheno.all.frq.modif','Filefreq1':'All/Freq/CIMT_All.pheno.nosmokers.frq.modif','Filefreq2':'All/Freq/CIMT_All.pheno.smokers.frq.modif','out':'All/all_imputed-meancIMT_Res.qassoc.final.gxe.withfreq'}
args = parseArguments()

if Test :
   FileGxE='All/all_imputed-meancIMT_Res.qassoc.final.gxe.modif'
   Filefreq1='All/Freq/CIMT_All.pheno.nosmokers.frq.modif'
   Filefreq2='All/Freq/CIMT_All.pheno.smokers.frq.modif'
   FilefreqAll = 'All/Freq/CIMT_All.pheno.all.frq.modif'
   out='Test'
else :
  FileGxE=args.file_gxe
  Filefreq1=args.freq1
  Filefreq2=args.freq2
  FilefreqAll=args.freqAll
  out=args.out

if args.dir_tmp==None:
   dirtmp="tmp/"
   try:
    # Create target Directory
     os.mkdir(dirtmp)
   except FileExistsError:
    print("Directory " , dirtmp,  " already exists")
else :
  dirtmp=args.dir_tmp

HeadI=dirtmp+"/"+os.path.basename(FileGxE)

gxe = pd.read_csv(FileGxE,delim_whitespace=True)
if args.bfile and args.pheno and args.pheno_file and args.pheno_gxe :
   datapheno=pd.read_csv(args.pheno_file, delim_whitespace=True) 
   gxe[args.gwas_rs].to_csv(HeadI+'list_rs.rs', header=False, index=False, mode='w',doublequote=False,sep='\t')
   datapheno=datapheno[datapheno[args.pheno].isna()==False] 
   lvalue=datapheno[args.pheno_gxe].unique().tolist() 
   lvalue=[x for x in lvalue if math.isnan(x)==False]
   if len(lvalue)!=2 :
      print(datapheno[args.pheno_gxe])
      print(lvalue)
      sys.exit('found different value unique in '+args.pheno_gxe)
   dataphenoshort=datapheno.iloc[:,[0,1]]
   Head=HeadI+"tmpall"
   dataphenoshort.to_csv(Head+'.ind', header=False, index=False, mode='w',doublequote=False,  sep='\t',)
   os.system("plink --keep-allele-order --threads "+args.plk_cores+" --keep "+Head+".ind -bfile "+args.bfile+" --freq --out "+Head+" --extract "+HeadI+"list_rs.rs")
   FilefreqAll=Head+'.frq'
   Head=HeadI+"temp_"+str(lvalue[0])
   datapheno1=datapheno[datapheno[args.pheno_gxe]==lvalue[0]]
   dataphenoshort=datapheno1.iloc[:,[0,1]]
   dataphenoshort.to_csv(Head+'.ind', header=False, index=False, mode='w',doublequote=False,  sep='\t',)
   os.system("plink --keep-allele-order --threads "+args.plk_cores+" --keep "+Head+".ind -bfile "+args.bfile+" --freq --out "+Head+" --extract "+HeadI+"list_rs.rs")
   Filefreq1=Head+'.frq'
   Head=HeadI+"temp_"+str(lvalue[1])
   datapheno1=datapheno[datapheno[args.pheno_gxe]==lvalue[1]]
   dataphenoshort=datapheno1.iloc[:,[0,1]]
   dataphenoshort.to_csv(Head+'.ind', header=False, index=False, mode='w',doublequote=False,  sep='\t',)
   os.system("plink --keep-allele-order --threads "+args.plk_cores+" --keep "+Head+".ind -bfile "+args.bfile+" --freq --out "+Head+" --extract "+HeadI+"list_rs.rs")
   Filefreq2=Head+'.frq'



   

freqAll = pd.read_csv(FilefreqAll,delim_whitespace=True)
freqAll = freqAll[['CHR','SNP', 'MAF','NCHROBS']]
freqAll['NCHROBS']=freqAll['NCHROBS']/2
freqAll['CHR']=freqAll['CHR'].replace(' ','').astype(str)
freqAll['SNP']=freqAll['SNP'].replace(' ','')
freqAll=freqAll.rename(index=str, columns={"MAF": "Freq_All", "NCHROBS": "N_All", "CHR":args.gwas_chr, "SNP":args.gwas_rs})
gxe[args.gwas_chr]=gxe[args.gwas_chr].astype(str)
gxeAll=pd.merge(gxe,freqAll,how="left",on=[args.gwas_chr, args.gwas_rs])
del gxe 
del freqAll
freq1 = pd.read_csv(Filefreq1,delim_whitespace=True)
freq1['NCHROBS']=freq1['NCHROBS']/2
freq1 = freq1[['CHR','SNP', 'MAF','NCHROBS']]
freq1['CHR']=freq1['CHR'].replace(' ','').astype(str)
freq1['SNP']=freq1['SNP'].replace(' ','')
freq1=freq1.rename(index=str, columns={"MAF": "Freq_1", "NCHROBS": "N_1", "CHR":args.gwas_chr, "SNP":args.gwas_rs})
gxeAll=pd.merge(gxeAll,freq1,how="left",on=[args.gwas_chr, args.gwas_rs])

del freq1

freq2 = pd.read_csv(Filefreq2,delim_whitespace=True)
freq2 = freq2[['CHR','SNP', 'MAF','NCHROBS']]
freq2['NCHROBS']=freq2['NCHROBS']/2
freq2['CHR']=freq2['CHR'].replace(' ','').astype(str)
freq2['SNP']=freq2['SNP'].replace(' ','')
freq2=freq2.rename(index=str, columns={"MAF": "Freq_2", "NCHROBS": "N_2", "CHR":args.gwas_chr, "SNP":args.gwas_rs})
gxeAll=pd.merge(gxeAll,freq2,how="left",on=[args.gwas_chr, args.gwas_rs])
del freq2

gxeAll.to_csv(out, sep='\t', na_rep='NA', header=True, index=False, mode='w')



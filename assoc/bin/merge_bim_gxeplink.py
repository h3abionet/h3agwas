#!/usr/bin/env python3
import re

# Merge PLINK GxE outpu with bim file 

import sys
import pandas as pd
import argparse
import numpy as np
import math


def GetInfo(readfile):
   a=readfile.readline().replace('\n','')
   return(a,a.split())

def parseArguments():
    parser = argparse.ArgumentParser(description='Merge PLINK GxE output with bim file')
    parser.add_argument('--plgxe',type=str,required=True, help="output of gxe of plink in same order than bim file")
    parser.add_argument('--bim',type=str,required=True,help="bim files")
    parser.add_argument('--out', type=str,help="output file",required=True)
    args = parser.parse_args()
    return args
args = parseArguments()

pos_rsbim=1
pos_posbim=3
pos_A1bim=4
pos_A2bim=5
pos_rsgxe=1

pos_Z=8
pos_se1=4
pos_se2=7

read_bim=open(args.bim)
read_gxe=open(args.plgxe)
write_out=open(args.out,'w')
write_out_notfind=open(args.out+".notfind",'w')
ent=re.sub('[ \t]+','\t', re.sub('^[ ]+','',read_gxe.readline().replace('\n','')))+"POS\tA1\tA2\tBetaGxE\tSeGxE\n"
write_out.write(ent)
write_out_notfind.write(ent)

(bimline,splbimline)=GetInfo(read_bim)
(gxeline,splgxeline)=GetInfo(read_gxe)


while bimline and gxeline :
   if splbimline[pos_rsbim]=='.':
      (bimline,splbimline)=GetInfo(read_bim)
      if not bimline :
        break
   if splgxeline[pos_Z]!="NA" and splgxeline[pos_se2]!="NA" and splgxeline[pos_se1]!="NA" :
        Se1=float(splgxeline[pos_se1])
        Se2=float(splgxeline[pos_se2])
        ZGxE=float(splgxeline[pos_Z])
        SeGxE=math.sqrt(Se1*Se1 + Se2*Se2)
        BetaGxE=ZGxE*SeGxE
   else :
        SeGxE="NA"
        BetaGxE="NA"
   if splgxeline[pos_rsgxe]=='.':
      write_out_notfind.write(re.sub('[ ]+','\t', re.sub('^[ ]+','',gxeline))+"\tNA\tNA\tNA"+str(BetaGxE)+"\t"+str(SeGxE)+"\n")
      (gxeline,splgxeline)=GetInfo(read_gxe)
      if not gxeline :
        break
   if splgxeline[pos_rsgxe]==splgxeline[pos_rsgxe]:
      write_out.write(re.sub('[ ]+','\t', re.sub('^[ ]+','',gxeline))+"\t"+splbimline[pos_posbim]+"\t"+splbimline[pos_A1bim]+"\t"+splbimline[pos_A2bim]+"\t"+str(BetaGxE)+"\t"+str(SeGxE)+"\n")
      (bimline,splbimline)=GetInfo(read_bim)
      (gxeline,splgxeline)=GetInfo(read_gxe)
   else :
      (bimline,splbimline)=GetInfo(read_bim)

for line in read_gxe :
    write_out_notfind.write(re.sub('[ ]+','\t', re.sub('^[ ]+','',line))+"\tNA\n")
 




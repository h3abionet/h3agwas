#!/usr/bin/env python3

# Merge PLINK GxE outpu with bim file 

import sys
import pandas as pd
import argparse
import numpy as np


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
pos_rsgxe=1
read_bim=open(args.bim)
read_gxe=open(args.plgxe)
write_out=open(args.out,'w')
write_out_notfind=open(args.out+".notfind",'w')
ent=read_gxe.readline().replace('\n','')+"\t"+"POS\n"
write_out.write(ent)
write_out_notfind.write(ent)

(bimline,splbimline)=GetInfo(read_bim)
(gxeline,splgxeline)=GetInfo(read_gxe)


while bimline and gxeline :
   if splbimline[pos_rsbim]=='.':
      (bimline,splbimline)=GetInfo(read_bim)
      if not bimline :
        break
   if splgxeline[pos_rsgxe]=='.':
      write_out_notfind.write(gxeline+"\tNA\n")
      (gxeline,splgxeline)=GetInfo(read_gxe)
      if not gxeline :
        break
   if splgxeline[pos_rsgxe]==splgxeline[pos_rsgxe]:
      write_out.write(gxeline+"\t"+splbimline[pos_posbim]+"\n")
      (bimline,splbimline)=GetInfo(read_bim)
      (gxeline,splgxeline)=GetInfo(read_gxe)
   else :
      (bimline,splbimline)=GetInfo(read_bim)

for line in read_gxe :
    write_out_notfind.write(line+"\tNA\n") 
 




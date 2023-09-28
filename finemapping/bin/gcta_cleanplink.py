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
TAB =chr(9)
#SNP A1 A2 freq b se p N 
def parseArguments():
    parser = argparse.ArgumentParser(description='format file for gcta, append N and frequencie if not present using bed file')
    parser.add_argument('--out', type=str,help="out of tex file",default="plink")
    parser.add_argument('--chr',type=str,help="n header in inp files",  required=True)
    parser.add_argument('--bin_plk',type=str,required=False,help="plink binary", default="plink")
    parser.add_argument('--bfile',type=str,required=True,help="bfile if need to compute frequency or N")
    parser.add_argument('--keep',type=str,required=False,help="file of data used for if need to compute frequency or N", default=None)
    parser.add_argument('--threads',type=int,required=False,help="", default=1)
    parser.add_argument('--pos_list',type=str,required=True,help="")
    args = parser.parse_args()
    return args

args = parseArguments()

fileposbed="a1234tmp.bed"
out=args.out
writebed=open(fileposbed,'w')
pos_list=[x for x in args.pos_list.split(',')]
chro=args.chr
for x in pos_list :
  writebed.write("\t".join([chro,x,x,chro+':'+x])+'\n')
writebed.close()
out=args.out
Cmd=args.bin_plk+" -bfile "+args.bfile+"  --keep-allele-order --make-bed "
if args.keep :
   Cmd+=" --keep "+args.keep
Cmd+=" --threads "+str(args.threads)
Cmd+=" --extract range  "+fileposbed
Cmd+=" --out "+ args.out
print(Cmd)
os.system(Cmd)

readbim=open(args.out+'.bim')
bimlist=[x.replace('\n','').split() for x in readbim]
readbim.close()
writebim=open(args.out+'.bim','w')
for x in bimlist :
 x[1]=x[0]+':'+x[3] 
 writebim.write("\t".join(x)+'\n')
writebim.close()

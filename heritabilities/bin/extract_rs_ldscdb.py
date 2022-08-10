#!/usr/bin/env python3

import sys
import pandas as pd
import argparse
import numpy as np
import gzip

EOL = chr(10)
#SNP A1 A2 freq b se p N 
def parseArguments():
    parser = argparse.ArgumentParser(description='format file for gcta, append N and frequencie if not present using bed file')
    parser.add_argument('--ldsc_input',type=str,required=True,help="", default=1)
    parser.add_argument('--gwas_input',type=str,required=True,help="", default=1)
    parser.add_argument('--out',type=str,required=True,help="", default=1)
    args = parser.parse_args()
    return args

args = parseArguments()

readposinput=open(args.gwas_input)
head=readposinput.readline()
dicchro={}
for line in readposinput :
  line=line.replace('\n','')
  spll=line.split()
  if spll[3] not in dicchro :
    dicchro[spll[3]]={}
  dicchro[spll[3]][spll[4]]=[line, 'NA']

readposinput.close()

readldsc=gzip.open(args.ldsc_input, 'rb')
head=readldsc.readline()
cmt=1
chro=-1
for line in readldsc :
  spl=line.decode('utf-8').split()
  if cmt>1 and chro != spl[0] :
     print('error differnt chro in ldsc file\n exit')
     sys.exit(-1)
     cmt=2
  chro=spl[0]
  if chro in dicchro :
    if spl[2] in  dicchro[chro] :
       dicchro[chro][spl[2]][1]=spl[1]

writeout=open(args.out, 'w')
for pos in dicchro[chro] : 
   if dicchro[chro][pos][1]!='NA' :
     writeout.write(dicchro[chro][pos][0]+'\t'+dicchro[chro][pos][1]+'\n') 
writeout.close()

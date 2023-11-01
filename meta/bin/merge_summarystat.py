#!/usr/bin/env python3

import sys
import os
import argparse
import math

def read_infofile(infofile) :
  readinfof=open(infofile)
  dicinfof={}
  for line in readinfof :
    spl=line.replace('\n','').split()
    dicinfof[spl[5]]=spl[0]+'\t'+spl[1]+'\t'+spl[2]+'\t'
  return dicinfof

def parseArguments():
    parser = argparse.ArgumentParser(description='transform file and header')
    parser.add_argument('--input_file',type=str,required=True, help="input file association")
    parser.add_argument('--out_file',type=str,required=True,help="output file")
    parser.add_argument('--info_file',type=str,required=True,help="list of header to print and replace and new header oldheader1:newheader1,oldheader2:newheader2")
    #MarkerName	Allele1	Allele2	Freq1	FreqSE	MinFreq	MaxFreq	Effect	StdErr	P-value	Direction	HetISq	HetChiSq	HetDf	HetPVal
    parser.add_argument('--metal',type=int,required=False, default=0)
    args = parser.parse_args()
    return args

args = parseArguments()

dicinfo=read_infofile(args.info_file)


readres=open(args.input_file)
writeres=open(args.out_file,'w')
header=readres.readline().replace('\n','')
if args.metal :
  splmetal=header.split()
  headmetal_dir=splmetal.index('Direction')
  writeres.write('SNP\tCHRO\tBP\tNStudies\t'+header+'\n')
else :
  writeres.write('SNP\tCHRO\tBP\t'+header+'\n')

for line in readres :
  if args.metal==1 :
    spl=line.replace('\n','').split()
    nbstudie=len(spl[headmetal_dir]) - spl[headmetal_dir].count("?")
    writeres.write(dicinfo[spl[0]]+str(nbstudie)+"\t"+line.replace('\n','')+'\n')
  else :
    writeres.write(dicinfo[spl[0]]+line.replace('\n','')+'\n')
  




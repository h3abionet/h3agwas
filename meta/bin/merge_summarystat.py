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
    args = parser.parse_args()
    return args

args = parseArguments()

dicinfo=read_infofile(args.info_file)


readres=open(args.input_file)
writeres=open(args.out_file,'w')
header=readres.readline().replace('\n','')
writeres.write('SNP\tCHRO\tBP\t'+header+'\n')
for line in readres :
  spl=line.split()
  writeres.write(dicinfo[spl[0]]+line.replace('\n','')+'\n')
  




#!/usr/bin/env python3

''' 
format file for gcta, append freq and n if need and bfile
'''
import argparse
import sys
import os

EOL = chr(10)
#SNP A1 A2 freq b se p N 
def parseArguments():
    parser = argparse.ArgumentParser(description='format file for gcta, append N and frequencie if not present using bed file')
    parser.add_argument('--inp_asso',type=str,required=True, help="association files")
    parser.add_argument('--out', type=str,help="out of tex file",default="test.tex")
    parser.add_argument('--chro_header',type=str,required=True,help="n header in inp files")
    parser.add_argument('--bp_header',type=str,help="bp header")
    args = parser.parse_args()
    return args


args = parseArguments()
inp     = args.inp_asso

bphead=args.bp_header
chrohead=args.chro_header

readfile=open(inp)

splline=readfile.readline().replace('\n','').split()
bpindex=splline.index(bphead)
chrindex=splline.index(chrohead)

writebed=open(args.out, 'w')
for line in readfile :
  spl=line.split()
  writebed.write(spl[chrindex]+"\t"+spl[bpindex]+"\t"+spl[bpindex]+"\t"+spl[chrindex]+':'+spl[bpindex]+'\n')




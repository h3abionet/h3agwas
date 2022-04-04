#!/usr/bin/env python3

# Converts PLINK covariate and fam file into a covariate file for Gemma

import sys
import argparse

def parseArguments():
    parser = argparse.ArgumentParser(description='fill in missing bim values')
    parser.add_argument('--bim',type=str,required=True)
    parser.add_argument('--out',type=str,required=True)
    args = parser.parse_args()
    return args

args=parseArguments()

readbim=open(args.bim)
writebim=open(args.out,'w')

for line in readbim :
  spll=line.replace('\n','').split() 
  a1=spll[4].upper()
  a2=spll[5].upper()
  newrs=spll[1]+'\t'+spll[0]+'_'+spll[3]
  if a1 > a2 :
    newrs+="_"+a1+"_"+a2
  else :
    newrs+="_"+a2+"_"+a1
  writebim.write(newrs+'\n')

      


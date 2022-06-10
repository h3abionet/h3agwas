#!/usr/bin/env python3

# Converts PLINK covariate and fam file into a covariate file for Gemma

import sys
import argparse
import re

def parseArguments():
    parser = argparse.ArgumentParser(description='fill in missing bim values')
    parser.add_argument('--bim',type=str,required=True)
    parser.add_argument('--out',type=str,required=True)
    parser.add_argument('--out_dup',type=str,required=True)
    args = parser.parse_args()
    return args

args=parseArguments()

readbim=open(args.bim)
writebim=open(args.out,'w')
writedup=open(args.out+'.dup','w')

listrs=set([])
for line in readbim :
  spll=line.replace('\n','').split() 
  a1=spll[4].upper()
  a2=spll[5].upper()
  a1=re.sub('[^A-Z]','',a1)
  a2=re.sub('[^A-Z]','',a2)
  if len(a1) > 5 :
     a1=a1[0:5]
  if len(a2) > 5 :
     a2=a2[0:5]
  newrs=spll[0]+'_'+spll[3]
  if a1 > a2 :
    newrs+="_"+a1+"_"+a2
  else :
    newrs+="_"+a2+"_"+a1
  cmtrs=0
  newrsi=newrs
  while newrs in listrs :
      writedup.write(spll[0]+'\t'+spll[3]+'\t'+spll[3]+'\t'+spll[0]+':'+spll[3])
      print("duplicate "+newrs)
      newrs=newrsi+'_'+str(cmtrs)
      cmtrs+=1
  listrs.add(newrs)
  writebim.write(spll[1]+'\t'+newrs+'\n')


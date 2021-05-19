#!/usr/bin/env python3

import sys
import pandas as pd
import argparse
import numpy as np

def readbed(filebed):
  readbed=open(filebed)
  dicbed={}
  for line in readbed:
     spll=line.replace('\n','').split()
     if spll[1] not in dicbed :
        dicbed[spll[1]]=[]
     dicbed[spll[1]].append([int(spll[2]),int(spll[3])])
  return dicbed

def checkrange(chro, pos, inforange):
     listblock=[]
     cmtblock=0
     if chro in inforange :
       for rangei in inforange[chro]:
        if pos>=rangei[0] and pos <=rangei[1]:
          listblock.append(cmtblock)
        cmtblock+=1
     return listblock


def parseArguments():
    parser = argparse.ArgumentParser(description='fill in missing bim values')
    parser.add_argument('--gwas',type=str,required=True)
    parser.add_argument('--bed',type=str,required=True)
    parser.add_argument('--chr_gwas', type=str,help="comma separated list of covariates",default="")
    parser.add_argument('--ps_gwas', type=str,help="comma separated list of covariates",default="")
    parser.add_argument('--out', type=str,help="output covariate file")
    args = parser.parse_args()
    return args


args=parseArguments()

bedres=readbed(args.bed)


readgwas=open(args.gwas)
headgw=readgwas.readline().replace('\n','').split()
gwas_bp=headgw.index(args.ps_gwas)
gwas_chr=headgw.index(args.chr_gwas)

writebloc=open(args.out, 'w')

for linegwas in readgwas :
   splgwas=linegwas.replace('\n','').split()
   chro=splgwas[gwas_chr]
   pos=int(splgwas[gwas_bp])
   listrange=checkrange(chro, pos, bedres)
   for block in listrange :
      #writebloc.write(chro+'_'+str(block)+'\t'+chro+'\t'+str(pos)+'\t'+str(bedres[chro][block][0])+'\t'+str(bedres[chro][block][1])+'\n')
      writebloc.write(chro+'_'+str(block)+'\t'+chro+'\t'+str(pos)+'\n')
     
writebloc.close()

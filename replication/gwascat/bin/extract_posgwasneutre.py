#!/usr/bin/env python3

# Converts PLINK covariate and fam file into a covariate file for Gemma

import sys
import argparse
def readbed(filebed, wind):
  readbed=open(filebed)
  dicpos={}
  dicrange={}
  for line in readbed:
     spll=line.replace('\n','').split()
     if spll[0] not in dicpos :
        dicpos[spll[0]]=[]
        dicrange[spll[0]]=[]
     bp=int(spll[1])
     dicpos[spll[0]].append(bp)
     dicrange[spll[0]].append([max(0,bp-wind),bp+wind,bp])
  return (dicpos, dicrange)

#    extract_posgwas.py --bed $pos --gwas $gwas --chr_gwas ${params.head_chr} --ps_gwas ${params.head_bp} --a1_gwas ${params.head_a1} --a2_gwas ${params.head_a2}  --wind ${params.size_win_kb}  --pval_gwas ${params.head_pval} --rs_gwas ${params.head_rs}
def parseArguments():
    parser = argparse.ArgumentParser(description='fill in missing bim values')
    parser.add_argument('--gwas',type=str,required=True)
    parser.add_argument('--bed',type=str,required=False)
    parser.add_argument('--chr_gwas', type=str,help="comma separated list of covariates",default="")
    parser.add_argument('--ps_gwas', type=str,help="comma separated list of covariates",default="")
    parser.add_argument('--wind', type=float,required=True,help="output covariate file")
    parser.add_argument('--pval_gwas', type=str,help="output covariate file")
    parser.add_argument('--out', type=str,help="output covariate file")
    args = parser.parse_args()
    return args


def checkpos(chro, pos, infopos):
     if chro in infopos and int(pos) in infopos[chro]:
        return True
     return False

def checkrange(chro, pos, inforange):
     if chro in inforange :
       for rangei in inforange[chro]:
        if pos>=rangei[0] and pos <rangei[1]:
         return True
     return False


args=parseArguments()

(infopos,inforange)=readbed(args.bed, (args.wind)*1000)
read_gwas=open(args.gwas)
gwashead=read_gwas.readline().replace('\n','')
gwasspl=gwashead.split()

chrgwas=gwasspl.index(args.chr_gwas)
posgwas=gwasspl.index(args.ps_gwas)
pvalgwas=gwasspl.index(args.pval_gwas)




writerange=open(args.out,'w')
writerange.write(gwashead+'\n')
for line in read_gwas :
  line=line.replace('\n','')
  spl=line.replace('\n','').split()
  chro=spl[chrgwas]
  pos=int(spl[posgwas])
  if checkrange(chro, pos, inforange)==False :
     writerange.write(line+'\n')

writerange.close()

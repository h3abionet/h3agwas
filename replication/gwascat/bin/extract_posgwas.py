#!/usr/bin/env python3

# Converts PLINK covariate and fam file into a covariate file for Gemma

import sys
import pandas as pd
import argparse
import numpy as np
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
    parser.add_argument('--a1_gwas',type=str,required=True,help="comma separated list of  pheno column")
    parser.add_argument('--a2_gwas', type=str,required=True,help="output fam file")
    parser.add_argument('--wind', type=float,required=True,help="output covariate file")
    parser.add_argument('--rs_gwas', required=True,type=str,help="output covariate file")
    parser.add_argument('--af_gwas', type=str,help="output covariate file")
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

(infopos,inforange)=readbed(args.bed, (args.wind+10)*1000)
read_gwas=open(args.gwas)
gwashead=read_gwas.readline().replace('\n','')
gwasspl=gwashead.split()
gwashead="\t".join(gwasspl)

chrgwas=gwasspl.index(args.chr_gwas)
posgwas=gwasspl.index(args.ps_gwas)
rsgwas=gwasspl.index(args.rs_gwas)
a1gwas=gwasspl.index(args.a1_gwas)
a2gwas=gwasspl.index(args.a2_gwas)
pvalgwas=gwasspl.index(args.pval_gwas)
#afgwas=gwasspl.index(args.a2_gwas)
if args.af_gwas :
   headmaf="FRQ"
   headmafpos=gwasspl.index(args.af_gwas)
   ListParam=[chrgwas,rsgwas, posgwas, a1gwas, a2gwas,headmafpos,pvalgwas, len(gwasspl)]
   ListHeadPlk=["CHR","SNP2", "BP", "A1", "A2", headmafchr,"P", "SNP"]
else :
   ListParam=[chrgwas,rsgwas, posgwas, a1gwas, a2gwas,pvalgwas, len(gwasspl)]
   ListHeadPlk=["CHR","SNP2", "BP", "A1", "A2", "P","SNP"]




writerange=open(args.out+'_range.init','w')
writerange_plk=open(args.out+'_range.assoc','w')
writerange_bed=open(args.out+'_range.bed','w')

writepos_bed=open(args.out+'_pos.bed','w')
writepos=open(args.out+'_pos.init','w')
writepos_plk=open(args.out+'_pos.assoc','w')
write_error=open(args.out+'_error.assoc','w')

#writeall_plk=open(args.out+'_all.assoc','w')

writepos.write(gwashead+'\tSNPplk'+'\n')
writerange.write(gwashead+'\tSNPplk'+'\n')


writepos_plk.write("\t".join(ListHeadPlk)+'\n')
writerange_plk.write("\t".join(ListHeadPlk)+'\n')

for line in read_gwas :
  line=line.replace('\n','')
  spl=line.replace('\n','').split()
  chro=spl[chrgwas]
  pos=int(spl[posgwas])
  if checkrange(chro, pos, inforange) :
     newrs=chro+"_"+str(pos)+"_"
     a1=spl[a1gwas].upper()
     a2=spl[a2gwas].upper()
     try :
       p=float(spl[pvalgwas])
     except :
       print(line)
       write_error.write('error pval not numeric '+spl[pvalgwas]+'\n')
       continue
     if a1 > a2 :
       newrs+=a1+"_"+a2
     else :
       newrs+=a2+"_"+a1
     spl.append(newrs)
     plkchar="\t".join([spl[x] for x in ListParam])
     newline="\t".join(spl)
     writerange.write(newline+'\n')
     writerange_plk.write(plkchar+'\n')
     writerange_bed.write(chro+"\t"+str(pos)+"\t"+str(pos)+'\t'+spl[rsgwas]+'\n')
     if checkpos(chro, pos, infopos):
         writepos.write(newline+'\n')
         writepos_plk.write(plkchar+'\n')
         writepos_bed.write(chro+"\t"+str(pos)+"\t"+str(pos)+'\t'+spl[rsgwas]+'\n')

writepos.close()
writerange.close()
writepos_plk.close()
writerange_plk.close()
#writeall_plk.close()

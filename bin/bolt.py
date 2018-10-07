#!/usr/bin/env python3
import os
import sys 
import re
cmd=" ".join(sys.argv[1::])+" 2> tmp.stderr"
print(sys.argv)
print(cmd)

os.system(cmd)
print(cmd)
lireerror=open("tmp.stderr")

BaliseEr=False
for lines in lireerror :
   if "ERROR: Heritability estimate" in lines :
      BaliseEr=True
   sys.stderr.write(lines)

cmdspl=re.split(r'[ =]', cmd)
ListStd=["--statsFile", "statsFileImpute2Snps"]
if BaliseEr:
   for CmtArg in range(len(cmdspl)):
       argv=cmdspl[CmtArg]
       for Std in ListStd :
          if Std in argv :
            File=cmdspl[CmtArg+1]
            Lire=open(File,'w')
            Lire.write("\t".join(["SNP","CHR","BP","GENPOS","ALLELE1","ALLELE0","A1FREQ","F_MISS","BETA","SE","P_BOLT_LMM_INF","P_BOLT_LMM"])+"\n")
            Lire.close

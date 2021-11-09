#!/usr/bin/env python3
import os
import time
import sys 
import re
FileStder="tmp.stderr"
Args=sys.argv
   
if "--reml" in " ".join(Args) :
    Cmt=0
    for x in sys.argv :
       if x == "--out_bolt2" :
         Balise=True 
         break
       Cmt+=1
    FileStdout=sys.argv[Cmt+1]
    Args.remove(FileStdout)
    Args.remove("--out_bolt2")
    cmd=" ".join(sys.argv[1::])+"> "+ FileStdout+" 2> "+ FileStder
          
else :
   cmd=" ".join(sys.argv[1::])+" 2> "+ FileStder
print(cmd)
error=os.system(cmd)
lireerror=open("tmp.stderr")

BaliseEr=False
AllLines=""
for lines in lireerror :
   if "ERROR: Heritability estimate" in lines :
      BaliseEr=True
   sys.stderr.write(lines)
   AllLines+=lines

cmdspl=re.split(r'[ =]', cmd)
ListStd=["statsFile", "statsFileImpute2Snps"]
  
if BaliseEr:
   if "--reml" in " ".join(sys.argv) :
      Ecrire=open(FileStdout,'w')
      Ecrire.write(AllLines)
      Ecrire.close()
   else :
     for CmtArg in range(len(cmdspl)):
       argv=cmdspl[CmtArg]
       for Std in ListStd :
          if Std in argv :
            File=cmdspl[CmtArg+1]
            WriteHead=open(File,'w')
            WriteHead.write("\t".join(["SNP","CHR","BP","GENPOS","ALLELE1","ALLELE0","A1FREQ","F_MISS","BETA","SE","P_BOLT_LMM_INF","P_BOLT_LMM"])+"\n")
            WriteHead.close()
else :
   for CmtArg in range(len(cmdspl)):
    argv=cmdspl[CmtArg]
    for Std in ListStd :
      if Std in argv :
         Read=open(cmdspl[CmtArg+1])
         Head=Read.readline()
         if len(Head)==0 :
           print("error file empty "+cmdspl[CmtArg+1])
           sys.exit(256) 
   if int(error)> 0:
     print(error)
     sys.exit(int(error))


#!/usr/bin/env python3


""" 

Takes input a relatedness file, a fam file, and a list of individuals and extracts the sub-matrix from the relatedness file
for the given individuals 


Jean-Tristan Brandenburg
"""


import sys
import pandas as pd
import numpy as np
import argparse

EOL=chr(10)

def errorMessage10(phe):
    print("""

    A problem has been detected in file <%s> column <%s>.

    There is some invalid data. I regret I can't tell you which row.


    Please check -- the data should be numeric only.


    If there is missing data, please use   NA



    """%(sys.argv[1],phe))

def parseArguments():
    parser = argparse.ArgumentParser(description='fill in missing bim values')
    parser.add_argument('--rel',type=str,required=True,help="File of relatdness matrix as gemma output")
    parser.add_argument('--phenofile',type=str,required=True,help="fam file use for compute relatdness matrix")
    parser.add_argument('--covfile',type=str,required=True,help="fam file use for compute relatdness matrix")
    parser.add_argument('--pospheno',type=int,required=True,help="fam file use for compute relatdness matrix")
    parser.add_argument('--relout',type=str,required=True,help="File with output pheno")
    parser.add_argument('--phenofileout',type=str,required=True,help="File with output pheno")
    parser.add_argument('--covfileout',type=str,required=True,help="File with output pheno")
    args = parser.parse_args()
    return args

args=parseArguments()

pospheno=args.pospheno
readpheno=open(args.phenofile)
NewHeader="FID IID\t"+readpheno.readline().split()[pospheno+1]
listeFIDKeep=[]
DicPheno={}
for Lines in readpheno :
    SplL=Lines.split()
    if SplL[1+pospheno]!='-9' and SplL[1+pospheno].upper()!="NA" :
       listeFIDKeep.append(SplL[0]+" "+SplL[1])
       DicPheno[SplL[0]+" "+SplL[1]]=SplL[0]+" "+SplL[1]+"\t"+SplL[1+pospheno]


readmat=open(args.rel)
linemat=readmat.readline()
listeFID=linemat.split('\t')
print(listeFID[1:5])
readmat.close()

ListePosKept=[0]
CmtFID=0
FinalIdList=[]
for FID in listeFID :
   if FID in listeFIDKeep :
      ListePosKept.append(CmtFID)
      FinalIdList.append(FID)
   CmtFID+=1 

readmat=open(args.rel)
writemat=open(args.relout, 'w')
CmtL=0
print('begin : open and write maatrix pheno in file '+args.relout)
for Line in readmat :
   Line=Line.replace('\n','')
   if CmtL in ListePosKept :
      Chaine=[]
      SplLine=Line.split('\t')
      for Pos in ListePosKept :
        Chaine.append(SplLine[Pos])
      writemat.write("\t".join(Chaine)+"\n")
   CmtL+=1
readmat.close()
writemat.close()
print('end : open and write maatrix pheno in file '+args.relout)

print('begin : write pheno in file '+args.phenofileout)
WritePheno=open(args.phenofileout,'w')
WritePheno.write(NewHeader+'\n')
for FID in FinalIdList :
    WritePheno.write(DicPheno[FID]+'\n')
WritePheno.close()
print('end : write pheno in file '+args.phenofileout)

readcov=open(args.covfile)
NewHeader=readcov.readline().replace('\n','')
DicCov={}
print('begin : red cov from file '+args.covfile)
for Lines in readcov :
    SplL=Lines.split()
    DicCov[SplL[0]+" "+SplL[1]]=Lines.replace('\n','')
readcov.close()
print('emd : red cov from file '+args.covfile)

print('begin : write cov '+args.covfileout)
writecov=open(args.covfileout, 'w')
writecov.write(NewHeader+'\n')
for FID in FinalIdList :
    writecov.write(DicCov[FID]+'\n')
writecov.close()
print('end : write cov')



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
    parser.add_argument('--inp_fam',type=str,required=True,help="fam file use for compute relatdness matrix")
    parser.add_argument('--relout',type=str,required=True,help="File with output pheno")
    parser.add_argument('--lind',type=str,required=True,help="liste ind to kept (one individual :FID IID by line)")
    args = parser.parse_args()
    return args

args=parseArguments()

readfam=open(args.inp_fam)
listeFID=[]
for Lines in readfam :
    SplL=Lines.split()
    listeFID.append(SplL[0]+"-"+SplL[1])
readfam.close()


readind=open(args.lind)
listeFIDKept=[]
for Lines in readind :
    SplL=Lines.split()
    listeFIDKept.append(SplL[0]+"-"+SplL[1])
readind.close()

ListePosKept=[]
CmtFID=0
for FID in listeFID :
   if FID in listeFIDKept :
      ListePosKept.append(CmtFID)
   CmtFID+=1 

readmat=open(args.rel)
writemat=open(args.relout, 'w')
CmtL=0
for Line in readmat :
   Chaine=[]
   if CmtL in ListePosKept :
      SplLine=Line.split()
      for Pos in ListePosKept :
        Chaine.append(SplLine[Pos])
      writemat.write("\t".join(Chaine)+"\n")
   CmtL+=1
readmat.close()
writemat.close()

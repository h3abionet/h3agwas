#!/usr/bin/env python3

import os
import argparse
import shutil
import sys

#      rename_mtag.py --listi $listf --listf $listf2 --dir_out mtag_res
def parseArguments():
    parser = argparse.ArgumentParser(description='format file for gcta, append N and frequencie if not present using bed file')
    parser.add_argument('--listi',type=str,required=True, help="association files")
    parser.add_argument('--listf', type=str,help="out of tex file",default="test.tex")
    parser.add_argument('--dir_out', type=str,help="out of tex file",default="test.tex")
    args = parser.parse_args()
    return args

args = parseArguments()
listI=args.listi.split(',')
listF=args.listf.split(',')
DirOut=args.dir_out
if len(listI)!=len(listF) :
   sys.exit("error len(listI)!=len(listF)")

##
for cmt in range(len(listI)) :
  headt="trait_"+str(cmt+1)+".txt" 
  FileMtag= DirOut+"/"+os.path.basename(listI[cmt])+".mtag"
  balise=False
  for File in listF:
     if headt in File :
        shutil.copy(File,FileMtag) 
        balise=True
  if balise==False:
     print(args.listi)
     print(args.listf)
     sys.exit('file with head '+headt+' not found\nexit') 


#!/usr/bin/env python3

import sys
import pandas as pd
import argparse
import numpy as np


EOL=chr(10)

def parseArguments():
    parser = argparse.ArgumentParser(description='fill in missing bim values')
    parser.add_argument('--list_file_annot',type=str,required=True, help="file contains chro and list of files with annotation")
    parser.add_argument('--info_pos',type=str,required=True,help="file contain one position : rs chro pos")
    parser.add_argument('--out', type=str,help="output")
    args = parser.parse_args()
    return args

args = parseArguments()

TAB =chr(9)

args = parseArguments()
read_infopos=open(args.info_pos)
infopos=read_infopos.readline().replace("\n","").split()
read_listanno=open(args.list_file_annot)
File=None
for x in read_listanno :
  tmpx=x.replace('\n','').split()
  if tmpx[0].upper()=="ALL" or tmpx[0]==infopos[1] :
    File=tmpx[1] 
    break

if File==None :
   print("doesn't find chromosome or All in "+ args.info_pos+"\n")
   sys.exit(1)

print("read : "+File+"\n")
read_info=open(File)
out=open(args.out, 'w')
out.write(read_info.readline())
for Line in read_info :
   Spl=Line.split()
   if infopos[1]==Spl[0] and infopos[2]==Spl[1] :
      out.write(Line)




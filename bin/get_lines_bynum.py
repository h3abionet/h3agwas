#!/usr/bin/env python3

import sys
import argparse


''' 
return lines in function of lines number 
prog --file file --lines a,b,c --out file
'''

def parseArguments():
    parser = argparse.ArgumentParser(description='fill in missing bim values')
    parser.add_argument('--file',type=str)
    parser.add_argument('--lines',type=str)
    parser.add_argument('--out',type=str)
    args = parser.parse_args()
    return args

args=parseArguments()
liste_numline=[int(x) for x in args.lines.split(",") if len(x)>0]
cmtline=0
read=open(args.file)
write=open(args.out, 'w')
for line in  read :
  if cmtline in liste_numline:
     write.write(line)
  cmtline+=1





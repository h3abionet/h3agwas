#!/usr/bin/env python3

import sys
import os
import argparse


def GetSep(Sep):
   if Sep==None :
     return None
   ListOfSep=["\t"," ",","]
   if len(Sep)>2 :
      Sep=Sep.upper()[:3]
      if Sep=='COM' :
         Sep=','
      elif Sep=='TAB' :
         Sep='\t'
      elif Sep=='WHI' :
         Sep=' '
   if Sep not in ListOfSep :
      return None
   return Sep



def parseArguments():
    parser = argparse.ArgumentParser(description='transform file and header')
    parser.add_argument('--input_file',type=str,required=True, help="input file association")
    parser.add_argument('--chro_header',type=str,required=True, help="input file association")
    parser.add_argument('--sep',type=str,default=None, required=False,help="if need to be at some chro")
    args = parser.parse_args()
    return args

args=parseArguments()
spl=GetSep(args.sep)
listchro=set([])
read=open(args.input_file)
posheadchro=read.readline().split(spl).index(args.chro_header)
for l in read :
   listchro.add(l.split(spl)[posheadchro])
print(" ".join(listchro))



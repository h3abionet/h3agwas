#!/usr/bin/env python3

import sys
import os
import argparse

def GetSep(Sep):
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
    parser.add_argument('--input_file',type=str,help="output file",required=True)
    parser.add_argument('--out_file',type=str,help="output file",required=True)
    parser.add_argument("--info_file", help="", type=str,required=True)
    args = parser.parse_args()
    return args

args=parseArguments()



infohead=args.info_file.split(",")
print (infohead)
l_infohead=[x.split(":")[0] for x in infohead]
l_filehead=[x.split(":")[1] for x in infohead]

nomrs=l_filehead[l_infohead.index('rsID')]
sep=GetSep(l_filehead[l_infohead.index('Sep')])

readfile=open(args.input_file)
writenew=open(args.out_file, 'w')

headeri=readfile.readline().replace('\n','').split(sep)
PosRs=headeri.index(nomrs)

writenew.write('rsID\n')
for ligne in readfile :
   spl=ligne.split(sep)
   if spl[PosRs]!='.' :
      writenew.write(spl[PosRs]+'\n')
writenew.close()

   





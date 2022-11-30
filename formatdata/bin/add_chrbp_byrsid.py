#!/usr/bin/env python3
import re
import sys
import os
import argparse
import gzip
import io

def isint(char) :
  try :
     int(char)
  except :
     return False
  return True

def read_chrbp(fileinfo) :
  dic={}
  readf=open(fileinfo)
  for line in readf :
      #9 205764 C T rs10811213 C
      spl=line.replace('\n','').split()
      dic[spl[2]]=[spl[0],spl[1]]
  return dic
#     add_chrbp_byrsid.py --gwas  $gwas --chrbp_info $outrs --rs_head ${params.head_rs} --chronew_head ${headnew_chr} --bpnew_head ${headnew_bp} --out $gwasf
def parseArguments():
    parser = argparse.ArgumentParser(description='transform file and header')
    parser.add_argument('--gwas',type=str,required=True, help="input file association")
    parser.add_argument('--chrbp_info',type=str,required=True, help="input file association")
    parser.add_argument('--rs_head',type=str,required=True,help="output file")
    parser.add_argument('--chronew_head',type=str,required=True,help="position of file to extract bp")
    parser.add_argument('--bpnew_head',type=str,required=True,help="rs of file to extract bp")
    parser.add_argument('--out',type=str, required=True)
    args = parser.parse_args()
    return args

args=parseArguments()
dicrs=read_chrbp(args.chrbp_info)
readgwas=open(args.gwas)
writegwas=open(args.out,'w')
writegwaserr=open(args.out+'.err','w')
head=readgwas.readline()
header=head.replace('\n','').split()
sep_out='\t'
rs_index=header.index(args.rs_head)
newhead=args.chronew_head+' '+args.bpnew_head+' '+head
writegwas.write(newhead)
writegwaserr.write(newhead)
for line in readgwas :
   spl=line.replace('\n','').split()
   rs=spl[rs_index] 
   chro=None
   bp=None
   if rs in dicrs :
     chro=dicrs[rs][0]
     bp=dicrs[rs][1]
   else :
     splrs=re.split("[:_]", rs)
     if len(splrs)>1 and isint(splrs[1]) :
       chro=splrs[0]
       bp=splrs[1]
   if chro and bp :
      writegwas.write(chro+sep_out+bp+sep_out+line) 
   else :
      writegwaserr.write(line) 


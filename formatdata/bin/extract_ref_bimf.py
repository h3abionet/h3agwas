#!/usr/bin/env python3

import sys
import os
import argparse
import gzip


def readfastagz(File) :
  Dic={}
  with gzip.open(File, "r") as ReadL :
   for ligne in ReadL:
    ligne=ligne.decode('utf-8').replace('\n','')
    if ligne[0]=='>' :
       Key=ligne.split(' ')[0].split('|')[0].split('\t')[0].replace('>','') 
       if Key in Dic :
         print(Dic.keys())
         print('fasta file '+ File+' more than one time same chro '+ Key)
         sys.exit(2)
       Dic[Key]=""
    Dic[Key]+=ligne
  return Dic

def readfastastdin() :
  Dic={}
  for ligne in sys.stdin:
    ligne=ligne.replace('\n','')
    if ligne[0]=='>' :
       Key=ligne.split(' ')[0].split('|')[0].split('\t')[0].replace('>','')
       if Key in Dic :
         print(Dic.keys())
         print('fasta file '+ File+' more than one time same chro '+ Key)
         sys.exit(2)
       Dic[Key]=""
    Dic[Key]+=ligne
  return Dic


def parseArguments():
    parser = argparse.ArgumentParser(description='transform file and header')
    parser.add_argument('--bim',type=str,required=True, help="input file association")
    parser.add_argument('--fasta',type=str,required=False, help="input file association")
    parser.add_argument('--out',type=str,required=True, help="input file association")
    parser.add_argument('--run', type=str, required=False, default='False')
    args = parser.parse_args()
    return args

args=parseArguments()
if args.run[0]=='T':
  if args.fasta :
    fasta=readfastagz(args.fasta)
  else :
    fasta=readfastastdin()
  readbim=open(args.bim)
  writetodel=open(args.out+'.del')
  writeref=open(args.out+'.allref')
  for line in readbim:
   spll=line.split()
   chro=spll[0]
   #0	200610-10	0	0	C	T
   #26	200610-37	0	16482	G	A
   if chro in fasta :
     pos=int(spll[3])
     all1=spll[4]
     all2=spll[5]
     allref=fasta[chro][pos-1]
     if allref == all1 or allref == all2:
         writeref.write(spll[1]+'\t'+allref+'\n') 
         writetodel.write("\t".join([chro,str(pos),spll[1],spll[4], spll[5], allref])+'\n')
   else :
    writetodel.write("\t".join([chro,str(pos),spll[1],spll[4], spll[5], "NA"])+'\n')

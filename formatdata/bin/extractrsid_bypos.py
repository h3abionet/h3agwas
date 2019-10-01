#!/usr/bin/env python3
import sys
import os
import argparse
import gzip
import io

def read_chrbp(file_chrbp) :
    read=open(file_chrbp)
    listchro=set([])
    for line in read :
      spl=line.split()
      listchro.add(spl[0]+" "+spl[1])  
    return listchro

def parseArguments():
    parser = argparse.ArgumentParser(description='transform file and header')
    parser.add_argument('--ref_file',type=str,required=True, help="input file association")
    parser.add_argument('--out_file',type=str,required=True,help="output file")
    parser.add_argument('--chro_ps',type=int,required=True,help="position of file to extract chro")
    parser.add_argument('--bp_ps',type=int,required=True,help="position of file to extract bp")
    parser.add_argument('--rs_ps',type=int,required=True,help="rs of file to extract bp")
    parser.add_argument('--chr',type=str,required=False,help="specific chromosome to extract")
    parser.add_argument('--file_chrbp',type=str, required=True)
    args = parser.parse_args()
    return args


args = parseArguments()

balisegz=False
if args.ref_file== 'stdin' :
  print("read in input file")
  read=sys.stdin 
else :
   #ext=os.path.splitext(args.ref_file)[1]
   #if ext in ['.gz','.gzip'] :
   #   print("read in gzip")
   #   read=gzip.open(args.ref_file,'rb')
   #   read=io.BufferedReader(read)
   #   balisegz=True
   #else :
   read=open(args.ref_file)

chrolist=read_chrbp(args.file_chrbp)
print("nb pos to found "+str(len(chrolist))+'\n')
write=open(args.out_file, 'w')
##
poschro=args.chro_ps
posbp=args.bp_ps
posrs=args.rs_ps

if args.chr :
  chro=args.chr
  for line in read :
     if line[0]=='#' :
        continue
     spl=line.split()
     if chro!=spl[poschro] :
        continue
     infopos=spl[poschro]+" "+spl[posbp]
     if infopos in chrolist :
       write.write(infopos+" "+spl[posrs]+'\n')
       chrolist.remove(infopos)
     if len(chrolist)==0 :
       break
else :
  for line in read :
     if line[0]=='#' :
        continue
     spl=line.split()
     infopos=spl[poschro]+" "+spl[posbp]
     if infopos in chrolist :
       write.write(infopos+" "+spl[posrs]+'\n')
       chrolist.remove(infopos)
     if len(chrolist)==0 :
       break

print("nb pos not found "+str(len(chrolist))+'\n')

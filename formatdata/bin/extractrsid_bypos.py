#!/usr/bin/env python3
import sys
import os
import argparse
import gzip
import io

def read_chrbp(file_chrbp) :
    read=open(file_chrbp)
    listchro=set([])
    listchro2={}
    for line in read :
      spl=line.split()
      listchro.add(spl[0]+" "+spl[1])  
      listchro2[spl[0]+" "+spl[1]]=[spl[2], spl[3]]
    return (listchro, listchro2)

def parseArguments():
    parser = argparse.ArgumentParser(description='transform file and header')
    parser.add_argument('--ref_file',type=str,required=True, help="input file association")
    parser.add_argument('--out_file',type=str,required=True,help="output file")
    parser.add_argument('--chro_ps',type=int,required=True,help="position of file to extract chro")
    parser.add_argument('--bp_ps',type=int,required=True,help="position of file to extract bp")
    parser.add_argument('--rs_ps',type=int,required=True,help="rs of file to extract bp")
    parser.add_argument('--a1_ps',type=int,required=True,help="a1 of file to extract bp")
    parser.add_argument('--a2_ps',type=int,required=True,help="a2 of file to extract bp")
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

(chrolist,chrodic)=read_chrbp(args.file_chrbp)
print("nb pos to found "+str(len(chrolist))+'\n')
write=open(args.out_file, 'w')
##
poschro=args.chro_ps
posbp=args.bp_ps
posrs=args.rs_ps
posa1=args.a1_ps
posa2=args.a2_ps

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
       spla1=spl[posa1].split(',')
       spla2=spl[posa2].split(',')
       if ((chrodic[infopos][0] in spla1) and (chrodic[infopos][1] in spla2)) or  ((chrodic[infopos][0] in spla2) and (chrodic[infopos][1] in spla1)) :
          write.write(infopos+" "+chrodic[infopos][0]+" "+chrodic[infopos][1]+" "+spl[posrs] +" "+spl[posa1]+'\n')
          chrolist.remove(infopos)
     #if len(chrolist)==0 :
     #  break
else :
  for line in read :
     if line[0]=='#' :
        continue
     spl=line.split()
     infopos=spl[poschro]+" "+spl[posbp]
     if infopos in chrolist :
       spla1=spl[posa1].split(',')
       spla2=spl[posa2].split(',')
       if (chrodic[infopos][0] in spla1 and chrodic[infopos][1] in spla2) or  (chrodic[infopos][0] in spla2 and chrodic[infopos][1] in spla1) :
         write.write(infopos+" "+chrodic[infopos][0]+" "+chrodic[infopos][1]+" "+spl[posrs]+" "+spl[posa1]+'\n')
         chrolist.remove(infopos)
     #if len(chrolist)==0 :
     #  break

print("nb pos not found "+str(len(chrolist))+'\n')

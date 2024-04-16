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
      spl=line.replace('\n','').split()
      listchro.add(spl[0]+" "+spl[1])  
      listchro2[spl[0]+" "+spl[1]]=[spl[2], spl[3], None]
    return (listchro, listchro2)

#3 3:61662:T:C 0 61662 C T
def  read_bim(bim) :
   definechr=""
   read=open(bim)
   listchro=set([])
   listchro2={}
   for line in read :
      spl=line.replace('\n','').split()
      listchro.add(spl[0]+" "+spl[3])
      if spl[1]=='.' :
         rs=None
      else :
         rs=spl[1]
      chro=spl[0]
      if "chr" in spl[0][0:3] :
         definechr="chr"
      if chro == '23' :
        chro=definechr+'X'
      listchro2[chro+" "+spl[3]]=[spl[4], spl[5], spl[1]]
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
    parser.add_argument('--file_chrbp',type=str, required=False, default=None)
    parser.add_argument('--bim',type=str, required=False, default=None)
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

if args.file_chrbp :
  (chrolist,chrodic)=read_chrbp(args.file_chrbp)
else :
  print('reading bim')
  (chrolist,chrodic)=read_bim(args.bim)

print("nb pos to found "+str(len(chrolist))+'\n')
write=open(args.out_file, 'w')
write2=open(args.out_file+'.init', 'w')
##
poschro=args.chro_ps
posbp=args.bp_ps
posrs=args.rs_ps
posa1=args.a1_ps
posa2=args.a2_ps
if args.bim :
  listrsnew=set([])
  listrsold=set([])
  write_rs=open(args.out_file+'.rs','w')

if args.chr :
  chro=args.chr
  for line in read :
     if line[0]=='#' :
        write2.write(line)
        continue
     spl=line.split()
     if chro!=spl[poschro] :
        continue
     infopos=spl[poschro]+" "+spl[posbp]
     if infopos in chrolist :
       spla1=spl[posa1].split(',')
       spla2=spl[posa2].split(',')
       if ((chrodic[infopos][0] in spla1) and (chrodic[infopos][1] in spla2)) or  ((chrodic[infopos][0] in spla2) and (chrodic[infopos][1] in spla1)) :
          write2.write(line)
          write.write(infopos+" "+chrodic[infopos][0]+" "+chrodic[infopos][1]+" "+spl[posrs] +" "+spl[posa1]+'\n')
          chrolist.remove(infopos)
          if args.bim and chrodic[infopos][2]:
             if (chrodic[infopos][2] not in listrsold) and (spl[posrs] not in listrsnew):
                write_rs.write(chrodic[infopos][2]+'\t'+spl[posrs]+'\t'+'\n')
                listrsold.add(chrodic[infopos][2])
                listrsnew.add(spl[posrs])
else :
  for line in read :
     if line[0]=='#' :
        write2.write(line)
        continue
     spl=line.split()
     infopos=spl[poschro]+" "+spl[posbp]
     if infopos in chrolist :
       spla1=spl[posa1].split(',')
       spla2=spl[posa2].split(',')
       if (chrodic[infopos][0] in spla1 and chrodic[infopos][1] in spla2) or  (chrodic[infopos][0] in spla2 and chrodic[infopos][1] in spla1) :
         write2.write(line)
         write.write(infopos+" "+chrodic[infopos][0]+" "+chrodic[infopos][1]+" "+spl[posrs]+" "+spl[posa1]+'\n')
         chrolist.remove(infopos)
         if args.bim and chrodic[infopos][2]:
             if (chrodic[infopos][2] not in listrsold) and (spl[posrs] not in listrsnew):
                write_rs.write(chrodic[infopos][2]+'\t'+spl[posrs]+'\t'+'\n')
                listrsold.add(chrodic[infopos][2])
                listrsnew.add(spl[posrs])

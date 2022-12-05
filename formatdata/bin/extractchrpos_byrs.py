#!/usr/bin/env python3
import sys
import os
import argparse
import gzip
import io
def GetSep(Sep):
   if not Sep :
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


def read_rs(filegwas, rs_head,a1_head, a2_head, splitgwas) :
    read=open(filegwas)
    head=read.readline().replace('\n', '').split(splitgwas)
    posrs=head.index(rs_head)
    posa1=head.index(a1_head)
    posa2=head.index(a2_head)
    listrs=set([])
    dicrs={}
    for line in read :
      spl=line.replace('\n','').split(splitgwas)
      listrs.add(spl[posrs]) 
      dicrs[spl[posrs]]=[spl[posa1].upper(),spl[posa2].upper()]
    return (listrs, dicrs)



def parseArguments():
    parser = argparse.ArgumentParser(description='transform file and header')
    parser.add_argument('--ref_file',type=str,required=True, help="input file association")
    parser.add_argument('--out_file',type=str,required=True,help="output file")
    parser.add_argument('--chro_ps',type=int,required=True,help="position of file to extract chro")
    parser.add_argument('--bp_ps',type=int,required=True,help="position of file to extract bp")
    parser.add_argument('--rs_ps',type=int,required=True,help="rs of file to extract bp")
    parser.add_argument('--a1_ps',type=int,required=True,help="a1 of file to extract bp")
    parser.add_argument('--a2_ps',type=int,required=True,help="a2 of file to extract bp")
    parser.add_argument('--gwas',type=str, required=True)
    parser.add_argument('--rs_head',type=str, required=True)
    parser.add_argument('--a1_head',type=str, required=True)
    parser.add_argument('--a2_head',type=str, required=True)
    parser.add_argument('--sep',type=str, required=False, default=None)
    args = parser.parse_args()
    return args


args = parseArguments()

balisegz=False
if args.ref_file== 'stdin' :
  print("read in input file")
  read=sys.stdin 
else :
   read=open(args.ref_file)

(list_rs,dicrs)=read_rs(args.gwas, args.rs_head, args.a1_head, args.a2_head,GetSep(args.sep))
print(list_rs)
write=open(args.out_file, 'w')
write2=open(args.out_file+'.init', 'w')
##
poschro=args.chro_ps
posbp=args.bp_ps
posrs=args.rs_ps
posa1=args.a1_ps
posa2=args.a2_ps

for line in read :
     if line[0]=='#' :
        write2.write(line)
        continue
     spl=line.split()
     infors=spl[posrs]
     if infors in list_rs:
       spla1=spl[posa1].split(',')
       spla2=spl[posa2].split(',')
       if (dicrs[infors][0] in spla1 and dicrs[infors][1] in spla2) or  (dicrs[infors][0] in spla2 and dicrs[infors][1] in spla1) :
         write2.write(line)
         write.write(spl[poschro]+" "+spl[posbp]+" "+spl[posrs]+'\n')


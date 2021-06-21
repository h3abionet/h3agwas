#!/usr/bin/env python3

import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import os

def ExtractLines(read_anno, write_anno, listpos) :
   CmtL=0
   lissave=set([])
   for Line in read_anno :
     Pos=int(Line.split()[1])
     if Pos in listpos and Line not in lissave:
        write_anno.write(Line) 
        lissave.add(Line)

def GetPosSearch(File) :
   read_listpos=open(File)
   listpos=set([])
   for Line in read_listpos :
     listpos.add(int(Line.split()[1]))
   read_listpos.close()
   return listpos

def GetPosSearch_int(File, write) :
   def search_pos(pos, listall):
      return [x[2] for x in listall if pos>=x[0] and pos<=x[1]]
   read_listpos=open(File)
   listpos=[]
   listline=[]
   for Line in read_listpos:
       pos=int(Line.split()[1])
       allinfo=search_pos(pos, listall)
       Chaine=Line.replace('\n', '')
       if len(allinfo) >0 :
          for cmt in range(len(allinfo[0])):
             listinfo=[]
             for cmt2 in range(len(allinfo)) :
                listinfo.append(allinfo[cmt2][cmt])
             Chaine+="\t"+";".join(list(set(listinfo)))
          write.write(Chaine+"\n")

def ReadIntervalFile(read_anno) :
   #listposbeg=[]
   #listposend=[]
   #listposinf=[]
   listall=[]
   for line in read_anno :
      spl=line.replace('\n','').split()
      listall.append([int(spl[1]),int(spl[2]),spl[3::]])
      #listposbeg.append(int(spl[1])) 
      #listposend.append(int(spl[2])) 
      #listposinfo.append(spl[3::])
   return sorted(listall)
       


def parseArguments():
    parser = argparse.ArgumentParser(description='Produces Manhatten, QQ plot and supporting tex file  for output from some tools')
    parser.add_argument('--list_pos',type=str,required=True, help="contains pos to extract (2 column)")
    parser.add_argument('--annov_file',type=str,required=True, help="annovar format file")
    parser.add_argument('--out', type=str,help="out of tex file",default="")
    args = parser.parse_args()
    return args

args = parseArguments()


#listpos.sort()
typeannov="f"
read_anno=open(args.annov_file)
tmphead=read_anno.readline().split()
balisef=True
if "Chr" not in tmphead[0] :
  typeannov="g"

Cmt=0
Cmt2=0
def isnum(x) :
  bal=True
  try :
    a=float(x)
  except :
    bal=False
  return bal

def gettype(x):
    if isnum(x) :
       xn=float(x) 
       if xn>=0 and xn<=1 :
         return 'freq'
       elif int(xn)==xn :
         if xn>1 and xn < 100000:
          return 'N'
         elif xn>0 :
          return 'bp'
         else :
          return 'on'
       else :
         return 'on'
    else :
      if x[0:2]=="rs" :
          return 'rs'
      else :
         x.upper()
         Ac=x.count('A')
         Tc=x.count('T')
         Cc=x.count('C')
         Gc=x.count('G')
         Nc=x.count('N')
         if Ac+Tc+Gc+Cc+Nc==len(x) :
           return 'all'
    return 'oc'

def getmax(ll):
   mval=max([ll[x] for x in ll])
   whichm=[x for x in ll if ll[x]==mval] 
   return whichm[0]

infotype=[{'rs':0, 'freq':0, 'N':0, 'bp':0, 'on':0, 'oc':0, 'all':0} for x in range(len(tmphead))]
for line in read_anno :
   spl=line.split()
   if isnum(spl[2]) and int(spl[2]) -int(spl[1])>1 :
     Cmt+=1
   Cmt2+=1
   if Cmt2>=1000:
      break
   for CmtCol in range(1,len(spl)):
        typebp=gettype(spl[CmtCol])
        infotype[CmtCol][typebp]+=1

print(infotype)
if Cmt/float(Cmt2) > 0.8 :
    typeannov="fx"

read_anno.close()
write_anno=open(args.out, 'w')
read_anno=open(args.annov_file)
if typeannov=="f" :
  print('f')
  listpos=GetPosSearch(args.list_pos)
  write_anno.write(read_anno.readline())
  ExtractLines(read_anno, write_anno, listpos)
elif typeannov=="g":
  print('g')
  listpos=GetPosSearch(args.list_pos)
  Ent=os.path.basename(args.annov_file).replace(".txt","")
  Tmp="#Chr\tStart"
  CmtNext=2
  if getmax(infotype[CmtNext])=='bp' :
     Tmp+="\tEnd"
     CmtNext+=1
  Tmp+="\tRef\tAlt"
  CmtNext+=2
  for cmcol in range(CmtNext, len(tmphead)) :
     print(getmax(infotype[cmcol]))
     Tmp+="\t"+getmax(infotype[cmcol])+"_"+Ent
  write_anno.write(Tmp+"\n")
  ExtractLines(read_anno, write_anno, listpos)
elif typeannov=="fx" :
    print('fx')
    head=read_anno.readline().replace('\n','').split()
    Tmp="#Chr\tStart\tEnd\tRef\tAlt\t"+"\t"+"\t".join(head[3::])+"\n"
    write_anno.write(Tmp)
    listall=ReadIntervalFile(read_anno)
    GetPosSearch_int(args.list_pos, write_anno)

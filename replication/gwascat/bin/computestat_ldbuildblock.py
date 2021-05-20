#!/usr/bin/env python3
''' 
manage output error for cojo in gcta
'''
from  subprocess import CalledProcessError
import subprocess
import os
import argparse
import numpy as np
import pandas as pd
import sys

EOL = chr(10)

def parseArguments():
    parser = argparse.ArgumentParser(description='extract from plink position ')
    parser.add_argument('--plink_ld',type=str,required=True, help="file error")
    parser.add_argument('--pos_cat',type=str,required=True, help="bed file, correponding to 3 column for input")
    parser.add_argument('--out',type=str,required=True, help="file error")
    args = parser.parse_args()
    return args

def readposcat(fileposcat):
  readpos=open(fileposcat) 
  dichropos={}
  for line in readpos :
      spl=line.replace('\n','').split()
      if spl[0] not in dichropos:
        dichropos[spl[0]]=set([])
      dichropos[spl[0]].add(int(spl[1]))  
  readpos.close()
  return dichropos
def readplinkld(plinkfile):
   print("begin :read ld from",plinkfile)
   readld=open(plinkfile)
   readld.readline()
   dicchro_ld={}
   for line in readld :
       splline=line.replace('\n','').split()
       chro=splline[0]
       if chro not in  dicchro_ld:
          print(" from ld chro", chro)
          dicchro_ld[chro]=[[],[],[]]
       pos1=int(splline[1])
       pos2=int(splline[4])
       dicchro_ld[chro][0].append(pos1)
       dicchro_ld[chro][1].append(pos2)
       dicchro_ld[chro][2].append(float(splline[6]))
   for chro in dicchro_ld.keys() :
       dicchro_ld[chro][0]= [x for _,x in sorted(zip(dicchro_ld[chro][2],dicchro_ld[chro][0]), reverse=True)] 
       dicchro_ld[chro][1]= [x for _,x in sorted(zip(dicchro_ld[chro][2],dicchro_ld[chro][1]), reverse=True)] 
       dicchro_ld[chro][2]=sorted(dicchro_ld[chro][2], reverse=True)
   print("end :read ld from",plinkfile)
   return dicchro_ld 
 
def getindice(ind, mylist):
  return [i for i, x in enumerate(mylist) if x == ind] 


args = parseArguments()
## keep alll position
dicposcat=readposcat(args.pos_cat)
infold=readplinkld(args.plink_ld)

cmtblok=0
dicbloc={}
writeld1=open(args.out+'_ld.out','w')
writeld2=open(args.out+'_ld_wind.out','w')

writeldmerg=open(args.out+'_ldext.out','w')
writeldmerg2=open(args.out+'_ldext_wind.out','w')

def checkwithotherbloc(bloc, listbloc):
 rlist=reversed(range(0, len(listbloc)))
 cmtblock2=0
 newlistadd=[]
 for bloc2 in listbloc:
    inters=bloc.intersection(bloc2)
    if len(inters)>0:
       newlistadd.append(cmtblock2)
    cmtblock2+=1
 return newlistadd
        

for chro in dicposcat.keys() :
   posnotfoundcat=[]
   dicbloc[chro]=[]
   if chro not in infold :
     continue 
   pos1ld=infold[chro][0]
   pos2ld=infold[chro][1]
   R2ld=infold[chro][2]
   listwind=[]
   cmtblock=0
   for poscat in dicposcat[chro]  :
     listposld=sorted(getindice(poscat,pos1ld)+getindice(poscat,pos2ld))
     if len(listposld)==0:
       posnotfoundcat.append([chro, poscat]) 
     else :
       bloc=set([])
       #dicbloc[chro]=set([])  
       for posld in listposld: 
         bloc.add(pos1ld[posld])
         bloc.add(pos2ld[posld])
       for posld in reversed(listposld) :
         pos1ld.pop(posld)
         pos2ld.pop(posld)
         R2ld.pop(posld)
       listpos=""
       for pos in bloc :
         writeld1.write(str(cmtblock)+'\t'+str(chro)+'\t'+str(pos)+'\n')
         listpos+=str(pos)+";"
       if len(bloc)>0 :
         minbloc=min(bloc)
         maxbloc=max(bloc)
         listwind.append([minbloc,maxbloc])
         writeld2.write(str(chro)+'_'+str(cmtblock)+'\t'+str(chro)+'\t'+str(minbloc)+'\t'+str(maxbloc)+'\t'+listpos+'\n')
       cmtblock+=1
       dicbloc[chro].append(bloc)
   for pos in posnotfoundcat :
     writeld1.write(str(cmtblock)+'\t'+str(chro)+'\t'+str(pos[1])+'\n')
     if len([x for x in listwind if pos[1]>=x[0] and pos[1]<=x[1]])==0 :
       writeld2.write(str(chro)+'_'+str(cmtblock)+'\t'+str(chro)+'\t'+str(pos[1])+'\t'+str(pos[1])+'\t'+str(pos[1])+'\n')
     cmtblock+=1
   balise=True 
   bloc=dicbloc[chro][0]
   dicbloc[chro].pop(0)
   newlistbloc=[]
   cmtblocn=1
   listwind=[]
   while balise :
     listbloctomerge=checkwithotherbloc(bloc, dicbloc[chro])
     if len(listbloctomerge)==0 :
        newlistbloc.append(bloc)
        bloc=dicbloc[chro][0]
        cmtblocn+=1
     else :
       listbloctomerge.sort(reverse=True)
       for posbloc in listbloctomerge :
         bloc.update(dicbloc[chro][posbloc])
         dicbloc[chro].pop(posbloc)
     if len(dicbloc[chro])==0:
       balise=False
   cmtblock=0
   for bloc in newlistbloc:
     listpos=""
     for pos in bloc :
       writeldmerg.write(str(chro)+')'+str(cmtblock)+'\t'+str(chro)+'\t'+str(pos)+'\n')
       listpos+=str(pos)+";"
     if len(bloc)>0 :
         minbloc=min(bloc)
         maxbloc=max(bloc)
         listwind.append([minbloc,maxbloc])
         writeldmerg2.write(chro+'_'+str(cmtblock)+'\t'+str(chro)+'\t'+str(minbloc)+'\t'+str(maxbloc)+'\t'+listpos+'\n')
     cmtblock+=1
   for pos in posnotfoundcat :
     writeldmerg.write(str(cmtblock)+'\t'+str(chro)+'\t'+str(pos[1])+'\n')
     if len([x for x in listwind if pos[1]>=x[0] and pos[1]<=x[1]])==0 :
       writeldmerg2.write(chro+'_'+str(cmtblock)+'\t'+str(chro)+'\t'+str(pos[1])+'\t'+str(pos[1])+'\t'+str(pos[1])+'\n')
     cmtblock+=1

writeld1.close()
writeld2.close()
     


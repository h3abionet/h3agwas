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
        dichropos[spl[0]]=[]
      dichropos[spl[0]].append(int(spl[1]))  
  readpos.close()
  return dichropos

def getblock(plink_ld) :
    def getblock(pos1, pos2,listbloc):
       cmtblock=0
       for block in listbloc :
         if (pos1 in block) or  (pos2 in block):
           return cmtblock
         cmtblock+=1
       return None
    def getblock(pos1, pos2,listbloc) :
       cmtblock=0
       listblock=[]
       for block in listbloc :
         if (pos1 in block) or  (pos2 in block):
           listblock.append(cmtblock)
         cmtblock+=1
       if len(listblock)==0 :   
         return None
       #if len(listblock)==1 :
       #  return listblock[0]
       #if len(listblock)>1 :
       #  print('multi block...')
       #  print(listblock)
       return listblock
 

    readplk=open(plink_ld)
    header=readplk.readline()
    chrocurrent=''
    dicchro={}
    for line in readplk : 
       splline=line.split() 
       chro=splline[0]
       if chro!=chrocurrent :
         if chrocurrent!='' :
           print('analyse of ld in ',chrocurrent,' finish. number of block',len(dicchro[chrocurrent]))
         chrocurrent=chro
       if chro not in dicchro:
         dicchro[chro]=[] 
       pos1=int(splline[1])
       pos2=int(splline[4])
       posblock=getblock(pos1, pos2,dicchro[chro])
       #print(pos1, pos2, posblock)
       if posblock==None :
          dicchro[chro]=[set([pos1,pos2])] + dicchro[chro]
       else :
         dicchro[chro][posblock[0]].add(pos1)
         dicchro[chro][posblock[0]].add(pos2)
         if len(posblock)>1 :
           newposblock=posblock[0]
           posblock.pop(0)
           posblock.sort()
           posblock.reverse()
           for cmb in posblock :
               dicchro[chro][newposblock].update(dicchro[chro][cmb])
               dicchro[chro].pop(cmb)
    readplk.close()
    return dicchro   
args = parseArguments()
## keep alll position
dicposcat=readposcat(args.pos_cat)
for chro in dicposcat.keys() :
    print("chro : "+chro+"pos number :"+str(len(dicposcat[chro])))
## group all ld in 
listld=getblock(args.plink_ld) 
## 

numblock=0
writeout=open(args.out,'w')
print( listld.keys())
allcat={}
for ldchro in listld.keys():
    allcat[ldchro]=[]
    if ldchro not in dicposcat:
       print('chro '+ldchro+' not in res of plink ld')
       listld.keys()
       continue
    for ldwind in listld[ldchro]:
       if len(ldwind)==1:
          print(ldwind)
          print("error block of 1")
       lisgwascat=[pos for pos in dicposcat[ldchro] if pos in ldwind]
       allcat[ldchro]+=lisgwascat
       if len(lisgwascat) > 0:
          print(len(ldwind),min(ldwind),max(ldwind))
          for pos in ldwind:
            if pos in lisgwascat:
               typepos="GC" 
            else :
               typepos="OT" 
            chaine="\t".join([str(numblock), str(ldchro),str(pos),typepos])+'\n'
            writeout.write(chaine)
          numblock+=1
           

for ldchro in dicposcat.keys() :
      typepos='CA'
      for pos in set(dicposcat[ldchro]) - set(allcat[ldchro])  :
        chaine="\t".join([str(numblock), str(ldchro),str(pos),typepos])+'\n'
        writeout.write(chaine)
        numblock+=1

writeout.close()

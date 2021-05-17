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
   print("read ld from",plinkfile)
   readld=open(plinkfile)
   readld.readline()
   dicchrold={}
   for line in readld :
       splline=line.replace('\n','').split()
       chro=splline[0]
       if chro not in  dicchrold:
          print(" from ld chro", chro)
          if chro!='1' :
           break
          dicchrold[chro]=[[],[],[]]
       pos1=int(splline[1])
       pos2=int(splline[4])
       dicchrold[chro][0].append(pos1)
       dicchrold[chro][1].append(pos2)
       dicchrold[chro][2].append(float(splline[6]))
   for chro in dicchrold.keys() :
       dicchrold[chro][0]= [x for _,x in sorted(zip(dicchrold[chro][2],dicchrold[chro][0]), reverse=True)] 
       dicchrold[chro][1]= [x for _,x in sorted(zip(dicchrold[chro][2],dicchrold[chro][1]), reverse=True)] 
       dicchrold[chro][2]=sorted(dicchrold[chro][2], reverse=True)
   return dicchrold 
 
def getindice(ind, mylist):
  return [i for i, x in enumerate(mylist) if x == ind] 


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
infold=readplinkld(args.plink_ld)

cmtblok=0
posnotfoundcat=[]
dicbloc={}
writeld1=open(args.out+'_ld.out','w')
writeld2=open(args.out+'_ld_resume.out','w')

writeldmerg=open(args.out+'_ldmerge.out','w')
writeldmerg2=open(args.out+'_ldmerge_resume.out','w')
cmtblock=0

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
   dicbloc[chro]=[]
   if chro not in infold :
     continue 
   pos1ld=infold[chro][0]
   pos2ld=infold[chro][1]
   R2ld=infold[chro][2]
   for poscat in dicposcat[chro]  :
     listposld=sorted(getindice(poscat,pos1ld)+getindice(poscat,pos2ld))
     if len(listposld)==0:
       posnotfoundcat.append([chro, poscat]) 
     else :
       print(poscat)
       bloc=set([])
       #dicbloc[chro]=set([])  
       for posld in listposld: 
         bloc.add(pos1ld[posld])
         bloc.add(pos2ld[posld])
       for posld in listposld :
         pos1ld.pop(posld)
         pos2ld.pop(posld)
         R2ld.pop(posld)
       listpos=""
       for pos in bloc :
         writeld1.write(str(cmtblock)+'\t'+str(chro)+'\t'+str(pos)+'\n')
         listpos+=str(pos)+";"
       if len(bloc)>0 :
         writeld2.write(str(cmtblock)+'\t'+str(chro)+'\t'+str(min(bloc))+'\t'+str(max(bloc))+'\t'+listpos+'\n')
       cmtblock+=1
       dicbloc[chro].append(bloc)
   balise=True 
   bloc=dicbloc[chro][0]
   dicbloc[chro].pop(0)
   newlistbloc=[]
   cmtblocn=1
   while balise :
     listbloctomerge=checkwithotherbloc(bloc, dicbloc[chro])
     print(listbloctomerge)
     if len(listbloctomerge)==0 :
        newlistbloc.append(bloc)
        bloc=dicbloc[chro][0]
        cmtblocn+=1
     else :
       listbloctomerge.sort(reverse=True)
       for posbloc in listbloctomerge :
         bloc.update(dicbloc[chro][posbloc])
         dicbloc[chro].pop(posbloc)
         print('bloc size running',len(bloc), len( newlistbloc))
     if len(dicbloc[chro])==0:
       balise=False
   cmtblock=0
   for bloc in newlistbloc:
     listpos=""
     for pos in bloc :
       writeldmerg.write(str(cmtblock)+'\t'+str(chro)+'\t'+str(pos)+'\n')
       listpos+=str(pos)+";"
     if len(bloc)>0 :
        writeldmerg2.write(str(cmtblock)+'\t'+str(chro)+'\t'+str(min(bloc))+'\t'+str(max(bloc))+'\t'+listpos+'\n')
     cmtblock+=1
     
writeld1.close()
writeld2.close()
     

#for chro in dicposcat.keys():
#   if chro not in infold :
#     continue 
#   for pos in dicposcat :
##for chro in dicposcat.keys() :
##    print("chro : "+chro+"pos number :"+str(len(dicposcat[chro])))
### group all ld in 
##listld=getblock(args.plink_ld) 
### 
 
#
##numblock=0
#writeout=open(args.out,'w')
#print( listld.keys())
#allcat={}
#for ldchro in listld.keys():
#    allcat[ldchro]=[]
#    if ldchro not in dicposcat:
#       print('chro '+ldchro+' not in res of plink ld')
#       listld.keys()
#       continue
#    for ldwind in listld[ldchro]:
#       if len(ldwind)==1:
#          print(ldwind)
#          print("error block of 1")
#       lisgwascat=[pos for pos in dicposcat[ldchro] if pos in ldwind]
#       allcat[ldchro]+=lisgwascat
#       if len(lisgwascat) > 0:
#          print(len(ldwind),min(ldwind),max(ldwind))
#          for pos in ldwind:
#            if pos in lisgwascat:
#               typepos="GC" 
#            else :
#               typepos="OT" 
#            chaine="\t".join([str(numblock), str(ldchro),str(pos),typepos])+'\n'
#            writeout.write(chaine)
#          numblock+=1
#           
#
#for ldchro in dicposcat.keys() :
#      typepos='CA'
#      for pos in set(dicposcat[ldchro]) - set(allcat[ldchro])  :
#        chaine="\t".join([str(numblock), str(ldchro),str(pos),typepos])+'\n'
#        writeout.write(chaine)
#        numblock+=1
#
#writeout.close()

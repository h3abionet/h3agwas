#!/usr/bin/env python3

import sys
import os
import argparse

def parseArguments():
    parser = argparse.ArgumentParser(description='transform file and header')
    parser.add_argument("--list_file", help="", type=str,required=True)
    parser.add_argument('--use_rs',type=str,help="if need to be limited at some rs", default=0)
    parser.add_argument("--out", help="output format ldsc, default none", type=str,required=True)
    args = parser.parse_args()
    return args

args=parseArguments()
splfile=args.list_file.split(',')
DicByRs={}
listRs=set([])
listChrBp={}
rsissue=''
listrsissue=set([])
listchrissue=set([])

for File in splfile :
   Fread=open(File)
   FreadL=Fread.readline().split()
   Fread.close()
   Fread=open(File)
   if len(FreadL)==3 :
     for line in Fread :
      splt=line.replace('\n', '').split()
      if splt[0] not in listRs :
         DicByRs[splt[0]]=[None,None,splt[1],splt[2],None]
      else :
         RsInfo=DicByRs[splt[0]]
         balisegood= (splt[1]==RsInfo[2] and splt[2]==RsInfo[3]) or  (splt[1]==RsInfo[3] and splt[2]==RsInfo[2])
         if balisegood ==False:
              listrsissue.add(splt[0]) 
      listRs.add(splt[0])            
   elif len(FreadL)==6:
       for line in Fread :
          splt=line.replace('\n', '').split() 
          NewRs=splt[5]
          if splt[0] not in listRs :  
               DicByRs[splt[0]]=[splt[1],splt[2],splt[3],splt[4], splt[5]]
          else :
             RsInfo=DicByRs[splt[0]]
             balisegood= (splt[3]==RsInfo[2] and splt[4]==RsInfo[3]) or  (splt[3]==RsInfo[3] and splt[4]==RsInfo[2])
             if balisegood ==False:
               listrsissue.add(splt[0]) 
             # check pos and chr
             if RsInfo[0] :
               if RsInfo[0] != splt[1] and RsInfo[1] != splt[2] :
                 listrsissue.add(splt[0]) 
             else :
                RsInfo[0]=splt[1] 
                RsInfo[1]=splt[2] 
                RsInfo[4]=splt[5] 
          listRs.add(splt[0])            
   else :
     print("colomn error number :"+str(len(FreadL)))
     sys.exit(3)
   print(len(listrsissue))
   print(len(listRs))


writeRs=open(args.out, 'w')
writeRs2=open(args.out+'_allinfo', 'w')
def checknone(x) :
  if x :
    return x
  return "NA"

for rs in DicByRs:
  RsInfo=DicByRs[rs]
  if rs not in listrsissue :
     if args.use_rs==1 : 
          writeRs.write(rs+'\t'+RsInfo[3]+'\t'+RsInfo[4]+'\n') 
     else :
       if RsInfo[0] :
         writeRs.write(rs+'\t'+'\t'.join([checknone(x) for x in RsInfo])+'\n') 
         writeRs2.write(rs+'\t'+'\t'.join([checknone(x) for x in RsInfo])+'\n') 
       else :
         listrsissue.add(rs)

writeRs.close()
writeRs2.close()

writeRsError=open(args.out+'_issue', 'w')
for rs in listrsissue :
    RsInfo=DicByRs[rs]
    writeRsError.write(rs+'\t'+'\t'.join([checknone(x) for x in RsInfo])+'\n')
writeRsError.close()

 

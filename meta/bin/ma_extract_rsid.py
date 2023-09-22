#!/usr/bin/env python3

import sys
import os
import argparse

def GetSep(Sep):
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





def parseArguments():
    parser = argparse.ArgumentParser(description='transform file and header')
    parser.add_argument('--input_file',type=str,help="output file",required=True)
    parser.add_argument('--out_file',type=str,help="output file",required=True)
    parser.add_argument("--info_file", help="", type=str,required=True)
    parser.add_argument("--ldsc", help="output format ldsc, default none",action="store_true")
    parser.add_argument("--out_sumstat",help="", type=str)
    args = parser.parse_args()
    return args

args=parseArguments()



infohead=args.info_file.split(",")
print (infohead)
l_infohead=[x.split(":")[0] for x in infohead]
l_filehead=[x.split(":")[1] for x in infohead]

nomrs=l_filehead[l_infohead.index('rsID')]
if "Sep" not in l_infohead :
   sep=None
else :
   sep=GetSep(l_filehead[l_infohead.index('Sep')])
if 'Chro' in l_infohead :
   nomchro=l_filehead[l_infohead.index('Chro')]
else :
   nomchro='NA'


if 'Pos' in l_infohead :
   nompos=l_filehead[l_infohead.index('Pos')]
else :
   nompos='NA'

A1Nom=l_filehead[l_infohead.index('A1')]
A2Nom=l_filehead[l_infohead.index('A2')]


readfile=open(args.input_file)
writenew=open(args.out_file, 'w')
writenewgwas=open(args.out_sumstat, 'w')

header_line=readfile.readline()
writenewgwas.write(header_line)
headeri=header_line.replace('\n','').split(sep)
PosRs=headeri.index(nomrs)
PosA1=headeri.index(A1Nom)
PosA2=headeri.index(A2Nom)
balisechro=False
if nompos.upper()!="NA" and nomchro.upper()!="NA" and len(nompos)>0 and len(nomchro)>0 :
   balisechro=True
   PosChro=headeri.index(nomchro)
   PosPos=headeri.index(nompos)
if args.ldsc :
   writenew2=open(args.out_file+'2', 'w')
   writenew.write('SNP\tA1\tA2\n')
   writenew2.write('SNP\tA1\tA2\tChro\tBP\n')
   if balisechro==False :
    for line in readfile :
     spl=line.replace('\n','').split(sep)
     spl[PosA1]=spl[PosA1].upper()
     spl[PosA2]=spl[PosA2].upper()
     if spl[PosRs]!='.' :
        writenew.write(spl[PosRs]+"\t"+spl[PosA1]+"\t"+spl[PosA2]+'\n')
   else :
    for line in readfile :
     spl=line.replace('\n','').split(sep)
     writenew.write(spl[PosRs]+"\t"+spl[PosA1]+"\t"+spl[PosA2]+'\n')
     pos=int(float(spl[PosPos]))
     if pos <=1 :
       print('warning : position <=0'+spl[PosPos]+'\n'+ '\t'.join(spl))
     else :
       spl[PosPos]=str(pos)
       writenew2.write(spl[PosRs]+"\t"+spl[PosA1]+"\t"+spl[PosA2]+'\t'+spl[PosChro]+'\t'+spl[PosPos]+'\n')
elif  balisechro == True  :
   writenew.write('rsID\tChro\tPos\tA1\tA2\tnewRs\n')
   for line in readfile :
     spl=line.replace('\n','').split(sep)
     spl[PosA1]=spl[PosA1].upper()
     spl[PosA2]=spl[PosA2].upper()
     if spl[PosA1] > spl[PosA2] :
         AA1=spl[PosA1]
         AA2=spl[PosA2]
     else :
         AA1=spl[PosA2]
         AA2=spl[PosA1]
     pos=int(float(spl[PosPos]))
     if pos <=1 :
       print('warning : position <=0'+spl[PosPos]+'\n'+ '\t'.join(spl))
     else :
       spl[PosPos]=str(pos)
       newrs=spl[PosChro]+'_'+spl[PosPos]+'_'+AA1+'_'+AA2
       if spl[PosRs]=='.' or spl[PosRs]=='NA' :
         spl[PosRs]=spl[PosChro]+':'+spl[PosPos] 
       writenewgwas.write('\t'.join(spl)+'\n')
       writenew.write(spl[PosRs]+"\t"+spl[PosChro]+"\t"+spl[PosPos]+"\t"+spl[PosA1]+"\t"+spl[PosA2]+"\t"+newrs+'\n')
else :
   writenew.write('rsID\tA1\tA2\n')
   for line in readfile :
     spl=line.replace('\n','').split(sep)
     if spl[PosRs]!='.' :
        writenew.write(spl[PosRs]+"\t"+spl[PosA1]+"\t"+spl[PosA2]+'\n')
        writenewgwas.write(line)
writenew.close()
writenewgwas.close()

   





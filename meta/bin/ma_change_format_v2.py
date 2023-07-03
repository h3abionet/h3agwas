#!/usr/bin/env python3

import sys
from scipy import stats
import os
import argparse
from math import sqrt, isnan

def formatrs(chro, bp,a1,a2) :
  a1=a1.upper()
  a2=a2.upper()
  if a1 > a2 :
     AA1=a1
     AA2=a2
  else :
     AA1=a2
     AA2=a1
  newrs=chro+'_'+bp+'_'+AA1+'_'+AA2
  return(newrs)

def reestimate_betase(pinit, betainit,n, freq,writemissing) :
    if pinit=='NA' or betainit=='NA' or n=='NA' or freq=='NA' :
       return ['NA', 'NA']
    if float(betainit)> 0 :
        sens=1
    else :
       sens = -1
    pinit=float(pinit)
    freq=float(freq)
    n = float(n)
    try :
      z=abs(stats.norm.ppf(1-pinit/2))*sens
      b=z/sqrt(2*freq*(1-freq)*(n+z**2))
      se=1/sqrt(2*freq*(1-freq)*(n+z**2))
      return [str(b),str(se)]
    except :
       writemissing.write('error computation '+ str(sens)+' '+ str(pinit) +' '+str(freq) + ' '+str(n)+' '+str(betainit)+'\n')
       return ["NA","NA"]


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

def ConfIsNotNull(x) :
   if x in ['','NA', 'Na', 'na'] :
      return False
   return True

def GetPosHead(head_file,head_tosearch):
   head_file=[x.replace(' ','')  for x in head_file]
   pos_head=[]
   for head in head_tosearch:
      if head not in head_file  :
         sys.exit("In "+ args.input_file+" not find "+  head+ "head \n exit")
      pos_head.append(head_file.index(head))
   return pos_head


def ChangeFile(File,FileOut ,old_header,new_header, sep, sepout):
    read=open(File)
    write=open(FileOut,'w')
    headinp=read.readline().replace('\n','').split(sep)
    pos_head=GetPosHead(headinp,old_header)
    write.write(sepout.join(new_header)+"\n")
    for line in read :
       spl=[x.replace(' ','') for x in line.replace('\n','').split(sep)]
       write.write(sepout.join([spl[x] for x in pos_head])+"\n")

def parseArguments():
    parser = argparse.ArgumentParser(description='transform file and header')
    parser.add_argument('--input_file',type=str,required=True, help="input file association")
    parser.add_argument('--out_file',type=str,required=True,help="output file")
    parser.add_argument('--info_file',type=str,required=True,help="list of header to print and replace and new header oldheader1:newheader1,oldheader2:newheader2")
    parser.add_argument('--sep_out',type=str,default="\t",help="separator output")
    parser.add_argument('--rs_ref',type=str,help="if need to be limited at some rs")
    parser.add_argument('--use_rs',type=str,help="if need to be limited at some rs", default=0)
    parser.add_argument('--used_pvalue',type=int,help="if need to be limited at some rs", default=0)
    args = parser.parse_args()
    return args
 


#ma_change_format.py  --input_file $file_assoc --out_file $newfile_assoc --info_file $info_file --rs_ref $file_ref
def checknull(x):
   if x=="" :
      return "NA"
   return x

## to defined if we need to add individual number
def GetCount(l_newhead,l_oldhead) :
  if ('Ncount' in l_newhead) :
    PosCount=l_newhead.index('Ncount')
    if   ('N' not in l_newhead):
        Ncount=l_oldhead[PosCount]
    else :
        Ncount=None
    del l_oldhead[PosCount]
    del l_newhead[PosCount]
  else :
    Ncount=None
  return (Ncount, l_newhead,l_oldhead)

args = parseArguments()

sep_out=GetSep(args.sep_out)


read=open(args.input_file)


infohead=args.info_file.split(",")
l_oldhead=[x.split(":")[1] for x in infohead if(ConfIsNotNull(x))]
l_newhead=[x.split(":")[0] for x in infohead  if(ConfIsNotNull(x))]

if 'Sep' in l_newhead :
   PosSep=l_newhead.index('Sep')
   sep=GetSep(l_oldhead[PosSep])
   del l_oldhead[PosSep]
   del l_newhead[PosSep]
else :
   sep =None

(Ncount,l_newhead,l_oldhead)=GetCount(l_newhead,l_oldhead)


rsId_inp=l_oldhead[l_newhead.index('rsID')]
balise_use_rs=False
balise_use_chrps=False
writemissing=open(args.out_file+'.error','w')
if args.rs_ref :
  open_rs=open(args.rs_ref)
  head=open_rs.readline()
  ncolrs=len(head.split())
  if ncolrs==1 :
     ls_rs=set([])
     for l in open_rs :
        rsspl=l.replace('\n', '').split()
        ls_rs.add(rsspl[0])
  else :
     ls_rs_dic={}
     ls_chrps_dic={}
     ls_rs=set([])
     for l in open_rs :
        rsspl=l.replace('\n', '').split()
        ls_rs.add(rsspl[0])
        ls_rs_dic[rsspl[0]]=[rsspl[1],rsspl[2], rsspl[3],rsspl[4], rsspl[5]]
        if rsspl[1] not in ls_chrps_dic :
          ls_chrps_dic[rsspl[1]]={} 
        try :
          ls_chrps_dic[rsspl[1]][int(rsspl[2])]=[rsspl,rsspl[5]]
        except :
           writemissing.write('error in position '+rsspl[1]+' '+rsspl[2] + ' '+rsspl+' in file '+args.rs_ref+'\n')
     if args.use_rs==1:
       ## will replace chro pos using rs
       balise_use_rs=True   
     else :
       ## will replace rs using chro pos
       balise_use_chrps=True

     
  
head_inp=read.readline().replace('\n','').split(sep)
ps_rsId_inp=head_inp.index(rsId_inp)
## sear
if balise_use_rs :
  if 'Pos' in l_newhead :
     HeadPosI=l_newhead.index('Pos')
     HeadPos=head_inp.index(l_oldhead[HeadPosI])
     del l_oldhead[HeadPosI]
     del l_newhead[HeadPosI]
  if 'Chro' in l_newhead :
     HeadChroI=l_newhead.index('Chro')
     HeadChro=head_inp.index(l_oldhead[HeadChroI])
     del l_oldhead[HeadChroI]
     del l_newhead[HeadChroI]

balise_use_dicrs=False
if balise_use_chrps :
  if 'Pos' in l_newhead :
     HeadPosI=l_newhead.index('Pos')
     HeadPos=head_inp.index(l_oldhead[HeadPosI])
  else :
    balise_use_dicrs=True
    balise_use_chrps=False
  if 'Chro' in l_newhead :
     HeadChroI=l_newhead.index('Chro')
     HeadChro=head_inp.index(l_oldhead[HeadChroI])
  else :
    balise_use_dicrs=True
    balise_use_chrps=False
  if 'rsID' in l_newhead :
    PosIndexI=l_newhead.index('rsID')
    del l_oldhead[PosIndexI]
    del l_newhead[PosIndexI]
    print("deleted index")

baliseN=False
if 'N' in l_newhead :
  baliseN=True

if Ncount:
  baliseN=True

balise_usedpvalue=False
if args.used_pvalue == 1 :
 balise_usedpvalue=True
 if (baliseN==False or ('freqA1' not in  l_newhead)  or ('Pval' not in l_newhead) or ('Beta' not in l_newhead)):
   print('baliseN', baliseN)
   print('freqA1 not in  l_newhead', 'freqA1' not in  l_newhead)
   print('Pval not in  l_newhead', 'Pval' not in  l_newhead)
   print('Beta not in  l_newhead', 'Beta' not in  l_newhead)
   print('used_pvalue must be have Beta, freq and value for N') 
   sys.exit(2)
 balise_usedpvalue=True
 if Ncount  :
   PosN=None
 else :
   PosN=head_inp.index(l_oldhead[l_newhead.index('N')])
 PosBetaI=l_newhead.index('Beta')
 PosBeta=head_inp.index(l_oldhead[PosBetaI])
 PosPval=head_inp.index(l_oldhead[l_newhead.index('Pval')])
 PosFreq=head_inp.index(l_oldhead[l_newhead.index('freqA1')])
 del l_oldhead[PosBetaI]
 del l_newhead[PosBetaI]
 if 'Se' in l_newhead:
    PosSeI=l_newhead.index('Se')
    del l_oldhead[PosSeI]
    del l_newhead[PosSeI]


print("balise_use_rs",balise_use_rs)
print("balise_use_chrps",balise_use_chrps)
print('balise_usedpvalue',balise_usedpvalue)

ps_head=GetPosHead(head_inp,l_oldhead)
print(ps_head, head_inp,l_oldhead)

balchangA1=False
balchangA2=False
if 'A1' in l_newhead :
  A1_inp=l_oldhead[l_newhead.index('A1')]
  ps_A1_inp=head_inp.index(A1_inp)
  balchangA1=True

if 'A2' in l_newhead :
  A2_inp=l_oldhead[l_newhead.index('A2')]
  ps_A2_inp=head_inp.index(A2_inp)
  balchangA2=True

listposfloat=[]
if  balise_usedpvalue  : 
  headplk=['RSID','CHR','BP','FREQ', 'P', 'N']
  headlist=['SNP','Chro','Pos','freqA1',  'Pval', 'N']
else :
  headplk=['RSID','CHR','BP','FREQ','BETA', 'SE', 'P', 'N']
  headlist=['SNP','Chro','Pos','freqA1', 'Beta', 'Se', 'Pval', 'N']

cmthead=0
l_newheadplk=[x for x in l_newhead]
for head in headlist:
   if head in l_newhead:
       indexlhead=l_newhead.index(head)
       l_newheadplk[indexlhead]=headplk[cmthead]
       newhead=l_oldhead[indexlhead]
       listposfloat.append(head_inp.index(newhead))
   cmthead+=1


l_newheadplk=[x.upper() for x in l_newheadplk]
if balise_use_rs and ('RSID' not in l_newheadplk):
   print(l_newheadplk)
   sys.exit('not rsid found')

if balise_use_chrps and ('CHR' not in l_newheadplk) and ('BP' not in l_newheadplk):
    sys.exit('chr and position column')

use_chrbp=False
if ('CHR' not in l_newheadplk) and ('BP' in l_newheadplk): 
  use_chrbp=True
  
p_minf=float('-inf')
p_pinf=float('inf')
def checkfloat(tmp, listposfloat):
   for x in listposfloat :
      try :
        resfl=float(tmp[x]) 
        if isnan(resfl) or resfl==p_minf or resfl==p_pinf:
           tmp[x]="NA"
      except ValueError:
        tmp[x]="NA"
   return tmp

write=open(args.out_file,'w')
writeplk=open(args.out_file+'.plk','w')
writeinfo=open(args.out_file+'.info','w')
## merge by rs
if balise_use_dicrs :
   print("balise_use_dicrs")
   l_newhead+=["rsID","CHRO", "POS"]
   l_newheadplk+=["SNP","CHR", "BP"]
   headtmp=[x.upper() for x in l_newhead]
   headtmpplk=[x.upper() for x in l_newheadplk]
   if Ncount :
      headtmp.append('N')
      headtmpplk.append('N')
   if balise_usedpvalue==1 :
      headtmp+=["BETA", "SE"]
      headtmpplk+=["BETA", "SE"]
   write.write(sep_out.join(headtmp)+"\n")
   writeplk.write(sep_out.join(headtmpplk)+"\n")
   for line in read :
     spl=line.replace('\n','').split(sep)
     if  spl[ps_rsId_inp] in ls_rs :
       spl=checkfloat(spl, listposfloat)
       if balchangA1 :
          spl[ps_A1_inp]=spl[ps_A1_inp].upper()
       if balchangA2 :
          spl[ps_A2_inp]=spl[ps_A2_inp].upper()
       tmp=[checknull(spl[x]) for x in ps_head]
       #ls_rs_dic[rsspl[0]]=[rsspl[1],rsspl[2], rsspl[3],rsspl[4]]
       tmprsdic=ls_rs_dic[spl[ps_rsId_inp]]
       tmp+=[tmprsdic[4],tmprsdic[0],tmprsdic[1]]
       if Ncount :
            Nval=Ncount
            tmp.append(Ncount)
       elif baliseN and balise_usedpvalue :
           Nval=spl[PosN]
       if balise_usedpvalue==1 :
           tmp+=reestimate_betase(spl[PosPval], spl[PosBeta],Nval, spl[PosFreq], writemissing)
       write.write(sep_out.join(tmp)+"\n")
       writeplk.write(sep_out.join(tmp)+"\n")
elif balise_use_rs :
   l_newhead+=["CHRO", "POS"]
   l_newheadplk+=["CHR", "BP"]
   headtmp=[x.upper() for x in l_newhead]
   headtmpplk=[x.upper() for x in l_newheadplk]
   if Ncount :
      headtmp.append('N')
      headtmpplk.append('N')
   if balise_usedpvalue==1 :
      headtmp+=["BETA", "SE"]
      headtmpplk+=["BETA", "SE"]
   write.write(sep_out.join(headtmp)+"\n")
   writeplk.write(sep_out.join(headtmpplk)+"\n")
   for line in read :
     spl=line.replace('\n','').split(sep)
     if  spl[ps_rsId_inp] in ls_rs :
       spl=checkfloat(spl, listposfloat)
       if balchangA1 :
          spl[ps_A1_inp]=spl[ps_A1_inp].upper()
       if balchangA2 :
          spl[ps_A2_inp]=spl[ps_A2_inp].upper()
       tmp=[checknull(spl[x]) for x in ps_head]
       tmp+=ls_rs_dic[spl[ps_rsId_inp]]
       if Ncount :
            Nval=Ncount
            tmp.append(Ncount)
       elif baliseN and balise_usedpvalue :
           Nval=spl[PosN]
       if balise_usedpvalue==1 :
           tmp+=reestimate_betase(spl[PosPval], spl[PosBeta],Nval, spl[PosFreq], writemissing)
       write.write(sep_out.join(tmp)+"\n")
       writeplk.write(sep_out.join(tmp)+"\n")
## balise use chro and positoin
elif balise_use_chrps :
   l_newhead+=['rsID']
   l_newheadplk+=['SNP']
   headtmp=[x.upper() for x in l_newhead]
   headtmpplk=[x.upper() for x in l_newheadplk]
   if Ncount :
      headtmp.append('N')
      headtmpplk.append('N')
   if balise_usedpvalue==1 :
      headtmp+=["BETA", "SE"]
      headtmpplk+=["BETA", "SE"]
   write.write(sep_out.join(headtmp)+"\n")
   writeplk.write(sep_out.join(headtmpplk)+"\n")
   for line in read :
     spl=line.replace('\n','').split(sep)
     bp=int(float(spl[HeadPos]))
     if  (spl[HeadChro] in ls_chrps_dic) and (bp in ls_chrps_dic[spl[HeadChro]]):
       desRS=ls_chrps_dic[spl[HeadChro]][bp]
       #spl[ps_rsId_inp]=ls_chrps_dic[spl[HeadChro]][spl[HeadPos]]
       spl=checkfloat(spl, listposfloat)
       if balchangA1 :
          spl[ps_A1_inp]=spl[ps_A1_inp].upper()
       if balchangA2 :
          spl[ps_A2_inp]=spl[ps_A2_inp].upper()
       tmp=[checknull(spl[x]) for x in ps_head]
       #tmp.append(ls_chrps_dic[spl[HeadChro]][spl[HeadPos]])
       newrs=formatrs(spl[HeadChro], str(bp), spl[ps_A1_inp], spl[ps_A2_inp])
       if newrs ==  desRS[1] :
         tmp.append(newrs)
         infors2="\t".join(desRS[0])+"\t"+newrs
         writeinfo.write(infors2+'\n')
         if Ncount :
            Nval=Ncount
            tmp.append(Ncount)
         elif baliseN and balise_usedpvalue :
           Nval=spl[PosN]
         if balise_usedpvalue==1 :
           tmp+=reestimate_betase(spl[PosPval], spl[PosBeta],Nval, spl[PosFreq], writemissing)
         write.write(sep_out.join(tmp)+"\n")
         writeplk.write(sep_out.join(tmp)+"\n")
       else :
           writemissing.write('error rs : '+newrs+' old rs '+desRS+'\n')
     else :
        writemissing.write('not found chr and bp :'+spl[HeadChro]+' '+str(bp)+'\n')
elif args.rs_ref :
   headtmp=[x.upper() for x in l_newhead]
   headtmpplk=[x.upper() for x in l_newheadplk]
   if Ncount :
      headtmp.append('N')
      headtmpplk.append('N')
   if balise_usedpvalue==1 :
      headtmp+=["BETA", "SE"]
      headtmpplk+=["BETA", "SE"]
   write.write(sep_out.join(headtmp)+"\n")
   writeplk.write(sep_out.join(headtmpplk)+"\n")
   for line in read :
     spl=line.replace('\n','').split(sep)
     spl=checkfloat(spl, listposfloat)
     if  spl[ps_rsId_inp] in ls_rs :
       if balchangA1 :
          spl[ps_A1_inp]=spl[ps_A1_inp].upper()
       if balchangA2 :
          spl[ps_A2_inp]=spl[ps_A2_inp].upper()
       tmp=[checknull(spl[x]) for x in ps_head]
       if Ncount :
          Nval=Ncount
          tmp.append(Ncount)
       elif baliseN and balise_usedpvalue :
          Nval=spl[PosN]
       if balise_usedpvalue==1 :
         tmp+=reestimate_betase(spl[PosPval], spl[PosBeta],Nval, spl[PosFreq], writemissing)
       write.write(sep_out.join(tmp)+"\n")
       writeplk.write(sep_out.join(tmp)+"\n")
else :
   tmphead=[x.upper() for x in l_newhead]
   tmpheadplk=[x.upper() for x in l_newheadplk]
   if Ncount :
      headtmp.append('N')
      headtmpplk.append('N')
   if balise_usedpvalue==1 :
      headtmp+=["BETA", "SE"]
      headtmpplk+=["BETA", "SE"]
   write.write(sep_out.join(tmphead)+"\n")
   writeplk.write(sep_out.join(tmpheadplk)+"\n")
   for line in read :
     spl=line.replace('\n','').split(sep)
     spl=checkfloat(spl, listposfloat)
     if balchangA1 :
          spl[ps_A1_inp]=spl[ps_A1_inp].upper()
     if balchangA2 :
          spl[ps_A2_inp]=spl[ps_A2_inp].upper()
     if use_chrbp :
        newrs=formatrs(spl[HeadChro], spl[HeadPos], spl[ps_A1_inp], spl[ps_A2_inp], writemissing)
     tmp=[checknull(spl[x]) for x in ps_head]
     if Ncount :
          Nval=Ncount
          tmp.append(Ncount)
     elif baliseN and balise_usedpvalue :
          Nval=spl[PosN]
     if balise_usedpvalue==1 :
       tmp+=reestimate_betase(spl[PosPval], spl[PosBeta],Nval, spl[PosFreq], writemissing)
     write.write(sep_out.join(tmp)+"\n")
     writeplk.write(sep_out.join(tmp)+"\n")

#!/usr/bin/env python3

import sys
import os
import argparse
import math

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
    args = parser.parse_args()
    return args
 


#ma_change_format.py  --input_file $file_assoc --out_file $newfile_assoc --info_file $info_file --rs_ref $file_ref
def checknull(x):
   if x=="" :
      return "NA"
   return x

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
rsId_inp=l_oldhead[l_newhead.index('rsID')]
baliserepchrpos=False
if(args.rs_ref) :
  open_rs=open(args.rs_ref)
  head=open_rs.readline()
  ncolrs=len(head.split())
  if ncolrs==1 :
     ls_rs=set([])
     for l in open_rs :
        rsspl=l.replace('\n', '').split()
        ls_rs.add(rsspl[0])
  else :
     baliserepchrpos=True
     ls_rs_dic={}
     ls_rs=set([])
     for l in open_rs :
        rsspl=l.replace('\n', '').split()
        ls_rs.add(rsspl[0])
        ls_rs_dic[rsspl[0]]=[rsspl[1],rsspl[2]]

## sear
if baliserepchrpos :
  if 'Pos' in l_newhead :
     PosPos=l_newhead.index('Pos')
     del l_oldhead[PosPos]
     del l_newhead[PosPos]
  
  if 'Chro' in l_newhead :
     PosPos=l_newhead.index('Chro')
     del l_oldhead[PosPos]
     del l_newhead[PosPos]

  
head_inp=read.readline().replace('\n','').split(sep)
ps_rsId_inp=head_inp.index(rsId_inp)
ps_head=GetPosHead(head_inp,l_oldhead)

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
headcmd=['freqA1', 'Beta', 'Se', 'Pval', 'N']
headplk=['FREQ','BETA', 'SE', 'P', 'N']

cmthead=0
for head in headcmd:
   if head in l_newhead:
       indexlhead=l_newhead.index(head)
       newhead=l_oldhead[indexlhead]
       l_newhead[indexlhead]=headplk[cmthead]
       listposfloat.append(head_inp.index(newhead))
   cmthead+=1
        

p_minf=float('-inf')
p_pinf=float('inf')
def checkfloat(tmp, listposfloat):
   for x in listposfloat :
      try :
        resfl=float(tmp[x]) 
        if math.isnan(resfl) or resfl==p_minf or resfl==p_pinf:
           tmp[x]="NA"
      except ValueError:
        tmp[x]="NA"
   return tmp

write=open(args.out_file,'w')
if baliserepchrpos :
   l_newhead+=["CHRO", "POS"]
   write.write(sep_out.join([x.upper() for x in l_newhead])+"\n")
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
       write.write(sep_out.join(tmp)+"\n")
elif args.rs_ref :
   write.write(sep_out.join([x.upper() for x in l_newhead])+"\n")
   for line in read :
     spl=line.replace('\n','').split(sep)
     spl=checkfloat(spl, listposfloat)
     if  spl[ps_rsId_inp] in ls_rs :
       if balchangA1 :
          spl[ps_A1_inp]=spl[ps_A1_inp].upper()
       if balchangA2 :
          spl[ps_A2_inp]=spl[ps_A2_inp].upper()
       write.write(sep_out.join([checknull(spl[x]) for x in ps_head])+"\n")
else :
   write.write(sep_out.join([x.upper() for x in l_newhead])+"\n")
   for line in read :
     spl=line.replace('\n','').split(sep)
     spl=checkfloat(spl, listposfloat)
     if balchangA1 :
          spl[ps_A1_inp]=spl[ps_A1_inp].upper()
     if balchangA2 :
          spl[ps_A2_inp]=spl[ps_A2_inp].upper()
     write.write(sep_out.join([checknull(spl[x]) for x in ps_head])+"\n")

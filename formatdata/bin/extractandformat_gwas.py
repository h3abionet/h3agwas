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
    write.close()
    read.close()

def parseArguments():
    parser = argparse.ArgumentParser(description='transform file and header')
    parser.add_argument('--input_file',type=str,required=True, help="input file association")
    parser.add_argument('--out_file',type=str,required=True,help="output file")
    parser.add_argument('--info_file',type=str,required=True,help="list of header to print and replace and new header oldheader1:newheader1,oldheader2:newheader2")
    parser.add_argument('--chr',type=str,default=None, required=False,help="if need to be at some chro")
    args = parser.parse_args()
    return args
 


#ma_change_format.py  --input_file $file_assoc --out_file $newfile_assoc --info_file $info_file --rs_ref $file_ref
def checknull(x):
   if x=="" :
      return "NA"
   return x

args = parseArguments()

writediscarded=open(args.out_file+'.discarded','w')


read=open(args.input_file)

infohead=args.info_file.split(",")
l_info=[x.split(':') for x in infohead]
l_infof=[]
l_oldheadf=[]
l_newheadf=[]

for x in l_info :
  if len(x)==3 and ConfIsNotNull(x[0]) and ConfIsNotNull(x[1]) and ConfIsNotNull(x[2]) :
    l_infof.append(x[0]) 
    l_oldheadf.append(x[1]) 
    l_newheadf.append(x[2]) 

if 'Sep' in l_infof :
   PosSep=l_infof.index('Sep')
   sep=GetSep(l_oldheadf[PosSep])
   newsep=GetSep(l_newheadf[PosSep])
   del l_oldheadf[PosSep]
   del l_infof[PosSep]
   del l_newheadf[PosSep]
   sep_out=newsep
else :
   sep =None
   sep_out='\t'


head_inp=read.readline().replace('\n','').split(sep)
nbcol_head=len(head_inp)
## sear
if 'Chro' in l_infof :
   pos_chro=head_inp.index(l_oldheadf[l_infof.index('Chro')])

balchangA1=False
if 'A1' in l_infof :
    balchangA1=True
    ps_A1_inp=head_inp.index(l_oldheadf[l_infof.index('A1')])

balchangA2=False
if 'A2' in l_infof :
    balchangA2=True
    ps_A2_inp=head_inp.index(l_oldheadf[l_infof.index('A2')])

ps_head=GetPosHead(head_inp,l_oldheadf)

listposfloat=[]
for head in ['freqA1', 'Beta', 'Se', 'Pval', 'N', 'Info']:
   if head in l_infof :
       newhead=l_oldheadf[l_infof.index(head)]
       listposfloat.append(head_inp.index(newhead))


p_minf=float('-inf')
p_pinf=float('inf')
def checkfloat(tmp, listposfloat):
   balise=True
   for x in listposfloat :
      try :
        resfl=float(tmp[x]) 
        if math.isnan(resfl) or resfl==p_minf or resfl==p_pinf:
           balise=False
           tmp[x]="NA"
      except ValueError:
        balise=False
        tmp[x]="NA"
   return (tmp, balise)

write=open(args.out_file,'w')
write.write(sep_out.join([x for x in l_newheadf])+"\n")
ChroSel=args.chr


if ChroSel :
  for line in read :
    spl=line.replace('\n','').split(sep)
    if spl[pos_chro]==ChroSel :
       if len(spl)!=nbcol_head :
         print('warning column number are not exact')
         writediscarded.write('NOTGODCol\t'+line)
         continue
       (spl,balgood)=checkfloat(spl, listposfloat)
       if balgood ==False:
          writediscarded.write("MISSING\t"+line)
          continue
       if balchangA1 :
          spl[ps_A1_inp]=spl[ps_A1_inp].upper()
       if balchangA2 :
          spl[ps_A2_inp]=spl[ps_A2_inp].upper()
       write.write(sep_out.join([checknull(spl[x]) for x in ps_head])+"\n")
else :
   for line in read :
     spl=line.replace('\n','').split(sep)
     if len(spl)!=nbcol_head :
         writediscarded.write('NOTGODCol\t'+line)
         print('warning column number are not exact')
         continue
     (spl,balgood)=checkfloat(spl, listposfloat)
     if balgood == False:
          writediscarded.write("MISSING\t"+line)
          continue
     if balchangA1 :
          spl[ps_A1_inp]=spl[ps_A1_inp].upper()
     if balchangA2 :
          spl[ps_A2_inp]=spl[ps_A2_inp].upper()
     write.write(sep_out.join([checknull(spl[x]) for x in ps_head])+"\n")

write.close()
writediscarded.close()

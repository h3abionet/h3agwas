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
write=open(args.out_file,'w')

infohead=args.info_file.split(",")
l_oldhead=[x.split(":")[1] for x in infohead if(ConfIsNotNull(x))]
l_newhead=[x.split(":")[0] for x in infohead  if(ConfIsNotNull(x))]

PosSep=l_newhead.index('Sep')
sep=GetSep(l_oldhead[PosSep])
del l_oldhead[PosSep]
del l_newhead[PosSep]
rsId_inp=l_oldhead[l_newhead.index('rsID')]
if(args.rs_ref) :
  open_rs=open(args.rs_ref)
  ls_rs=set([])
  for l in open_rs :
     ls_rs.add(l.replace('\n', '').split()[0])

head_inp=read.readline().replace('\n','').split(sep)
ps_rsId_inp=head_inp.index(rsId_inp)
## sear


ps_head=GetPosHead(head_inp,l_oldhead)
write.write(sep_out.join([x.upper() for x in l_newhead])+"\n")
balise=True
for line in read :
  spl=line.replace('\n','').split(sep)
  if args.rs_ref : 
     if spl[ps_rsId_inp] in ls_rs :
        balise=True
     else :
        balise=False
  if balise :
     write.write(sep_out.join([checknull(spl[x]) for x in ps_head])+"\n")


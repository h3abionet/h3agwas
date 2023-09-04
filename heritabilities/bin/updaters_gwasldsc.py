#!/usr/bin/env python3

import sys
import pandas as pd
import argparse
import numpy as np
import gzip

EOL = chr(10) 
#SNP A1 A2 freq b se p N 
def parseArguments():
    parser = argparse.ArgumentParser(description='format file for gcta, append N and frequencie if not present using bed file')
    parser.add_argument('--gwas',type=str,required=True,help="", default=1)
    parser.add_argument('--rstoupdate',type=str,required=True,help="", default=1)
    parser.add_argument('--rs_header',type=str,required=True,help="rs header in inp files")
    parser.add_argument('--chro_header',type=str,help="n header in inp files", default="")
    parser.add_argument('--bp_header',type=str,required=True,help="beta header in inp files")
    parser.add_argument('--p_header',type=str,required=True,help="p-value header in inp files")
    parser.add_argument('--n_header',type=str,required=False,help="p-value header in inp files")
    parser.add_argument('--info_header',type=str,required=False,help="p-value header in inp files")
    parser.add_argument('--a1_header',type=str,required=True,help="a1 header in inp files")
    parser.add_argument('--a2_header',type=str,required=True,help="a2 header in inp files")
    parser.add_argument('--af_header',type=str,required=True,help="a2 header in inp files")
    parser.add_argument('--beta_header',type=str,required=True,help="a2 header in inp files")
    parser.add_argument('--se_header',type=str,required=True,help="a2 header in inp files")
    parser.add_argument('--out',type=str,required=True,help="", default=1)
    parser.add_argument('--out_pos',type=str,required=True,help="", default=1)
    args = parser.parse_args()
    return args
args=parseArguments()

def gethead(header, listhead) :
  listposhead=[]
  for head in listhead :
   if head not in header : 
      print('error : '+head+" not in header of "+args.gwas+" ")
      print("list of head "+ " ".join(header))
      sys.exit(-1)
   listposhead.append(header.index(head))
  return listposhead


def getkey(chro, bp, a1, a2) :
   key=chro+str(bp)
   if a1 > a2 :
      key=key+'_'+a1+"_"+a2
   else :
      key=key+'_'+a2+"_"+a1
   return key

if args.out == args.out_pos :
  print("--out == --out_pos")
  sys.exit(-1)
       
read_rsupd=open(args.rstoupdate)
dic_rs={}
for line in  read_rsupd :
  #9:205764:C:T    T       C       9       205764  NA
  (oldrs, a1,a2,chro,bp, newrs)=line.replace('\n','').split() 
  if newrs !='NA' :
     key= getkey(chro, bp, a1, a2)
     dic_rs[key]=newrs 

read_rsupd.close()

listkeyrs=set(dic_rs)
readgwas=open(args.gwas)
headgwas=readgwas.readline()
splhead=headgwas.replace('\n','').split()

(rsidx, a1idx, a2idx, chridx, bpidx, rsidx, betaidx, seidx, pidx ,afidx)=gethead(splhead, [args.rs_header, args.a1_header, args.a2_header,args.chro_header, args.bp_header, args.rs_header, args.beta_header, args.se_header, args.p_header,args.af_header])
balisen=False
if args.n_header:
  nidx=gethead(splhead, [args.n_header])
  nidx=nidx[0]
  balisen=True

baliseinfo=False
if args.info_header :
  ninfox=gethead(splhead, [args.info_header])
  ninfox=ninfox[0]
  print(ninfox)
  baliseinfo=True
 

writegwas=open(args.out, 'w')
head=['SNP','CHR', 'BP','A1','A2','BETA', 'SE','P','FRQ']

if balisen :
   head.append('N')  

if baliseinfo :
   head.append('INFO')

writegwas.write(" ".join(head)+"\n")
writegwas2=open(args.out_pos, 'w')
writegwas2.write('SNP\tA1\tA2\n')
for line in  readgwas :
  spll=line.replace('\n','').split()
  spll[a1idx]=spll[a1idx].upper()
  spll[a2idx]=spll[a2idx].upper()
  newkey=getkey(spll[chridx], spll[bpidx], spll[a1idx], spll[a2idx])
  if newkey in listkeyrs :
    spll[rsidx]=dic_rs[newkey] 
    linegwas=[spll[rsidx],spll[chridx], spll[bpidx], spll[a1idx], spll[a2idx],spll[betaidx],spll[seidx],spll[pidx],spll[afidx]]
    if balisen :
      linegwas.append(spll[nidx])
    if baliseinfo :
      linegwas.append(spll[ninfox])
    writegwas.write(" ".join(linegwas)+'\n')
    writegwas2.write(spll[rsidx]+"\t"+spll[a1idx]+"\t"+spll[a2idx]+"\n")

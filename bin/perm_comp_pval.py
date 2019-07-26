#!/usr/bin/env python3

# computed p-value with permuted data with two way by locus and with minimum by simulation 

TAB=chr(9)
import sys
import pandas as pd
import argparse
import numpy as np

#SNP A1 A2 freq b se p N 
def parseArguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--inp',type=str,required=True,help="a2 header in inp files")
    parser.add_argument('--listgwas',type=str,required=True,help="se header in inp files")
    parser.add_argument('--out',type=str,help="n header in inp files", default=None)
    parser.add_argument('--head_pv',type=str,required=True,help="header for pvalue")
    #parser.add_argument('--head_rs',type=str,required=True,help="header for rs")
    args = parser.parse_args()
    return args

#https://groups.google.com/forum/#!topic/plink2-users/hT6qQs2b1VI
def GetPvalue(listfile,listpvlocalobs, headpv) :
    listminpv=[]#*len(listpvlocal) 
    listpvloc=[0]*len(listpvlocalobs)
    listpvlocalobs=[float(x) for x in listpvlocalobs]
    CmtF=0
    for File in listfile :
       Read=open(File)
       spllin=Read.readline().split()
       PosPva=spllin.index(headpv)	
       listpvloca=[float(line.split()[PosPva]) for line in Read] 
       listminpv.append(min(listpvloca))
       for x in range(len(listpvlocalobs)) :
           if listpvloca[x]<=listpvlocalobs[x] :
               listpvloc[x]+=1
       if CmtF%100==1:
          print(CmtF)
       CmtF+=1
    listminpval=[]
    listpvloc=[(x+1)/float(len(listminpv)+1) for x in listpvloc]
    for x in listpvlocalobs :
        listminpval.append((len([y for y in listminpv if y <= x])+1)/float(len(listminpv)+1))
    return (listpvloc, listminpval, listminpv)
args = parseArguments()

read = open(args.listgwas)
listfile=[x.replace("\n","") for x in read.readlines()]
read.close()


datagwas=pd.read_csv(args.inp,delim_whitespace=True)
(listpvloc, listminpva, listminpv)=GetPvalue(listfile, datagwas[args.head_pv].tolist(), args.head_pv)
datagwas['PvalmaxT']=listminpva
datagwas['PvalByLoc']=listpvloc
datagwas.to_csv(args.out,sep=TAB,header=True,index=False,na_rep="NA")
writemin=open(args.out+".MINPV", 'w')
writemin.write("\n".join([str(x) for x in listminpv]))



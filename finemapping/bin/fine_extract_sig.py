#!/usr/bin/env python3

# Converts PLINK covariate and fam file into a covariate file for Gemma

import sys
import pandas as pd
import argparse
import numpy as np
import os
import scipy.stats as stats


def extractinfobim(chro, begin, end,bimfile):
   readbim=open(bimfile) 
   listpos=[]
   listref=[]
   listalt=[]
   listsnp=[]
   listdup=set([])
   for line in readbim :
      splline=line.split()
      pos=int(splline[3])
      if chro == splline[0] and pos>=begin and pos<=end and (pos not in listdup):
         if pos in listpos :
           posindex=listpos.index(pos)
           listpos.pop(posindex) 
           listref.pop(posindex) 
           listalt.pop(posindex) 
           listdup.add(pos)
           listsnp.pop(posindex)
         else : 
          listpos.append(pos)
          listsnp.append(splline[1])
          listref.append(splline[4])
          listalt.append(splline[5])
   return (listpos, listref, listalt, listsnp)

def appendfreq(bfile, result, freq_header,rs_header, n_header, chr_header,bp_header, a1_header, a2_header, bin_plk, keep, threads) :
   plkfreqfil=os.path.basename(bfile)
   out_range="tmp.frq"
   listcol=[chr_header,bp_header,bp_header,rs_header]
   small[listcol].to_csv(out_range, sep=TAB, header=False, index=False,na_rep="NA")
   if bfile==None :
     print("no header for n or freq and bfile")
     sys.exit()
   Cmd=bin_plk+" -bfile "+bfile+" --freq --keep-allele-order --extract  range  "+out_range
   if keep :
     Cmd+=" --keep "+keep
   Cmd+=" --threads "+str(threads)
   Cmd+=" --out "+plkfreqfil
   os.system(Cmd)
   data_n=pd.read_csv(plkfreqfil+".frq",delim_whitespace=True)
   data_n.columns= ["CHR","SNP","A1_tmp","A2_tmp","MAF","NCHROBS"]
   databim = pd.read_csv(bfile+".bim", delim_whitespace=True)
   #CHR            SNP   A1   A2
   databim.columns = ["CHR","SNP", "Other","BP","A1_tmp", "A2_tmp"]
   data_n=data_n.merge(databim,how="inner", on=["CHR","SNP", "A1_tmp", "A2_tmp"])
   # CHR            SNP   A1   A2          MAF  NCHROBS
   #SNP A1 A2 freq b se p N
   data_n['N']=data_n['NCHROBS']/2
   if n_header==None and freq_header==None:
      data_n=data_n[["CHR","BP","MAF","N","A1_tmp",'A2_tmp']]
      freq_head="MAF"
      n_header="N"
   elif n_header==None :
      data_n=data_n[["CHR","BP","N"]]
      n_header="N"
   elif freq_header==None :
      data_n=data_n[["CHR","BP","MAF",'A1_tmp','A2_tmp']]
      freq_head="MAF"
   IsFreq=True
   IsN=True
   result[chr_header]=result[chr_header].astype(str)
   data_n['CHR']=data_n['CHR'].astype(str)
   result=result.merge(data_n,how="inner", left_on=[chr_header,bp_header], right_on=["CHR","BP"])
   balise=result[a2_header] == result['A1_tmp']
   result.loc[balise,'MAF']= 1 - result.loc[balise,'MAF']
   return (freq_head, n_header, result)

EOL=chr(10)

def parseArguments():
    parser = argparse.ArgumentParser(description='extract rs in gwas file')
    parser.add_argument('--inp_resgwas',type=str,required=True)
    parser.add_argument('--rs',type=str,required=False, help="list rs", default="")
    parser.add_argument('--chro',type=str,required=False, help="list rs", default="")
    parser.add_argument('--begin',type=float,required=False, help="list rs")
    parser.add_argument('--end',type=float,required=False, help="list rs")
    parser.add_argument('--min_pval',type=float,required=False, help="list rs")
    parser.add_argument('--chro_header',type=str,required=True,help="chro header in inp files")
    parser.add_argument('--n_header',type=str,required=False,help="chro header in inp files")
    parser.add_argument('--pos_header',type=str,required=True,help="pos header in inp files")
    parser.add_argument('--beta_header',type=str,required=True,help="beta header in inp files")
    parser.add_argument('--se_header',type=str,required=True,help="beta header in inp files")
    parser.add_argument('--p_header',type=str,required=True,help="beta header in inp files")
    parser.add_argument('--rs_header',type=str,required=False,help="beta header in inp files")
    parser.add_argument('--a1_header',type=str,required=True,help="beta header in inp files")
    parser.add_argument('--around_rs',type=float,required=False,help="beta header in inp files")
    parser.add_argument('--a2_header',type=str,required=True,help="beta header in inp files")
    parser.add_argument('--freq_header',type=str,required=False,help="frequencies header in inp files")
    parser.add_argument('--maf',type=float,default=0.0,help="minor allele frequencies")
    parser.add_argument('--out_head',type=str,default="out",help="around rs (pb)")
    parser.add_argument('--bin_plk',type=str,required=False,help="plink binary", default="plink")
    parser.add_argument('--bfile',type=str,required=False,help="bfile if need to compute frequency or N", default=None)
    parser.add_argument('--keep',type=str,required=False,help="file of data used for if need to compute frequency or N", default=None)
    parser.add_argument('--z_pval',type=int,required=False,help="file of data used for if need to compute frequency or N", default=0)
    parser.add_argument('--threads',type=int,required=False,help="", default=1)
    parser.add_argument('--n',required=False, help="bim file ")
    args = parser.parse_args()
    return args


args = parseArguments()

result = pd.read_csv(args.inp_resgwas,delim_whitespace=True, dtype={args.chro_header:str})
result[args.chro_header]=result[args.chro_header].astype(str)
result.drop_duplicates(subset=[args.chro_header, args.pos_header],keep=False,inplace=True)
rs=args.rs
freq_header=args.freq_header
n_header=args.n_header
rs_header=args.rs_header
maf=args.maf
out_file=args.out_head+"_finemap.z"
out_fileZ=args.out_head+"_caviar.z"
if args.rs :
   sub_result=result[result[args.rs_header]==rs]
   chro=sub_result[args.chro_header].tolist()[0]
   pos=sub_result[args.pos_header].tolist()[0]
   begin=pos-args.around_rs
   end=pos+args.around_rs
else :
   chro=args.chro
   begin=args.begin
   end=args.end
TAB=chr(9)
##rsid chromosome position allele1 allele2 maf beta se


(listbim, listalt, listref, listsnp)=extractinfobim(chro, begin,end,args.bfile+".bim")

small=result[(result[args.chro_header]==chro) & (result[args.pos_header]>=begin) & (result[args.pos_header]<=end)]

small=pd.merge(small, pd.DataFrame({"pos":listbim, 'altbim232':listalt, 'refbim232':listref, 'snpbim232':listsnp}),left_on=args.pos_header, right_on="pos" )
small=small[((small['altbim232']==small[args.a1_header]) & (small['refbim232']==small[args.a2_header])) | ((small['altbim232']==small[args.a2_header]) & (small['refbim232']==small[args.a1_header]))]

if args.min_pval and small[args.p_header].min()> args.min_pval:
   sys.exit('min pval '+str(args.min_pval)+'>'+' min(p) '+ str(small[args.p_header].min()))
small=small.loc[small[args.pos_header].isin(listbim)]
if args.n : 
  small['N']=args.n
  n_header='N'

if (n_header==None or freq_header==None):
  (freq_header, n_header, small)=appendfreq(args.bfile, small, freq_header,rs_header, n_header, args.chro_header,args.pos_header,args.a1_header, args.a2_header,args.bin_plk, args.keep, args.threads)
print(small)

PosCol=['snpbim232',args.chro_header, args.pos_header, args.a1_header, args.a2_header, freq_header, args.beta_header, args.se_header, args.p_header, n_header]
RsName=["rsid","chromosome","position","allele1","allele2","maf", "beta", "se", "p","N"]
small=small[(small[rs_header]==rs) | ((small[freq_header]>maf) & (small[freq_header]<(1-maf)))]

small=small[PosCol] 
small.columns=RsName
small=small.sort_values(["position"])
### for gcta
smallgcta=small[["rsid","chromosome","position","allele1","allele2","maf", "beta", "se", 'p','N']]
smallgcta=smallgcta.rename(columns={"rsid": "SNP", "chromosome": "chr", "position":'bp', 'allele1':'A1', 'allele2':'A2', 'beta':'b', 'maf':'freq'})

out_gcta=args.out_head+'.gcta'
smallgcta[['SNP','A1','A2','freq','b','se','p','N']].to_csv(out_gcta, sep=TAB, header=True, index=False,na_rep="NA")
maxn=smallgcta['N'].max()
writen=open('n.out', 'w')
writen.write(str(maxn))
writen.close()

## change freq and allele
bal=small['maf']>0.5
small['allele1_tmp']=small['allele1']
small['allele2_tmp']=small['allele2']
small['maf_tmp']=small['maf']
small['beta_tmp']=small['beta']
small.loc[bal,'allele1']=small.loc[bal,'allele2_tmp']
small.loc[bal,'allele2']=small.loc[bal,'allele1_tmp']
small.loc[bal,'beta']= - small.loc[bal,'beta_tmp']
small.loc[bal,'maf']= 1 - small.loc[bal,'maf_tmp']

if args.z_pval==0 :
  small['Z']=small['beta']/small['se']
else :
  tmpbeta=small['beta'].copy()
  tmpbeta[tmpbeta>=0]=1
  tmpbeta[tmpbeta<=0]=-1
  small['Z']=stats.norm.ppf(1-small['p']/2)
  small['Z']=small['Z'].abs()*tmpbeta
  small['beta']=small['Z']*small['se']

#small.to_csv('test_z', sep=' ', header=True, index=False,na_rep="NA")

small[["rsid","chromosome","position","allele1","allele2","maf", "beta", "se"]].to_csv(out_file, sep=" ", header=True, index=False,na_rep="NA")
small[["rsid","Z"]].to_csv(out_fileZ, sep=" ", header=False, index=False,na_rep="NA")

out_range=args.out_head+".range"
small[["chromosome","position","position", 'rsid']].to_csv(out_range, sep=TAB, header=False, index=False,na_rep="NA")

out_all=args.out_head+".all"
small.to_csv(out_all, sep=TAB, header=True, index=False,na_rep="NA")

out_range=args.out_head+".paintor"
small[["Z"]].to_csv(out_range, sep=TAB, header=True, index=False,na_rep="NA")
out_pos=args.out_head+'.pos'
small[["chromosome","position"]].to_csv(out_pos, sep=" ", header=True, index=False,na_rep="NA")

out_pos=args.out_head+'.rs'
small[["rsid"]].to_csv(out_pos, sep=" ", header=True, index=False,na_rep="NA")



#!/usr/bin/env python3
import argparse
from math import sqrt
import os
import sys
from scipy import stats
EOL=chr(10)
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

def extractinfobim(bimfile):
   readbim=open(bimfile)
   dicrskey={}
   dickeyrs={}
   listkey=set([])
   listchrbp=set([])
   listdup=set([])
   for line in readbim :
      splline=line.replace('\n','').split()
      chro=splline[0]
      rs=splline[1]
      bp=splline[3]
      a1=splline[4]
      a2=splline[5]
      keyposalt=formatrs(chro, bp,a1,a2)
      keychr=chro+':'+bp
      if keychr in listchrbp :
         listdup.add(keychr)
      else :
         listchrbp.add(keychr)
      listkey.add(keyposalt) 
      dicrskey[rs]=keyposalt
      dickeyrs[keyposalt]=rs
   return (listkey, listdup, dicrskey,dickeyrs)

def readsummarystat(sumstat, dickeyrs,listdup, listkey, rs_header,n_header, chr_header, bpval_header, beta_header, z_header, a1_header, a2_header, se_header, af_header, pval_header,n_value, maf, nlim, keep_genet) :
    def getposheader( head, headersumstat):
        nheader=None
        if head is not None:
          head=head.lower()
          try :
           nheader=headersumstat.index(head)
          except :
            print(headersumstat)
            print(head+" not found in "+" ".join(headersumstat))
            sys.exit(2)
        print('head '+str(head)+' found in column '+str(nheader))
        return nheader
    def gv(infogwas, posheader, othervalue="NA", formatfct=str) :
       if posheader is not None:
          try :
             val=infogwas[posheader]
          except :
            print(infogwas, posheader, othervalue, formatfct)
            sys.exit(2)
          try :
            return formatfct(val)
          except :
            return othervalue
       return othervalue
    readsummstat=open(sumstat)
    headersumstat=[x.lower() for x in readsummstat.readline().replace('\n','').split()]
    n_headerp=getposheader(n_header, headersumstat)
    se_headerp=getposheader(se_header, headersumstat)
    a2_headerp=getposheader(a2_header, headersumstat)
    a1_headerp=getposheader(a1_header, headersumstat)
    z_headerp=getposheader(z_header, headersumstat)
    beta_headerp=getposheader(beta_header, headersumstat)
    bpval_headerp=getposheader(bpval_header, headersumstat)
    chr_headerp=getposheader(chr_header, headersumstat)
    rs_headerp=getposheader(rs_header, headersumstat)
    af_headerp=getposheader(af_header, headersumstat)
    pval_headerp=getposheader(pval_header, headersumstat)
    if maf is not None:
      mafupper= 1 - maf
    if n_value is None:
       n_value_default="NA"
    dicres={}
    listnewkey=set([])
    balisen=(nlim  is not None) and (n_headerp is not None)
    ncol=len(headersumstat)
    for line in readsummstat :
        splline=line.replace('\n','').split()
        if len(splline)!=ncol :
           print(splline)
           print("column number different between header and line\n")
           continue
        key=formatrs(splline[chr_headerp],splline[bpval_headerp],splline[a1_headerp],splline[a2_headerp])
        key2=splline[chr_headerp]+':'+splline[bpval_headerp]
        balisersbim=key in dickeyrs
        if (key2 not in listdup)  and (keep_genet==0 or balisersbim):
           freq=gv(splline,af_headerp)
           balise=True
           if freq!="NA" and (maf is not None):
             freq=float(freq)
             balise=freq>=maf and freq<=mafupper
           if balisen and balise:
             Nval=float(splline[n_headerp])
             balise=balise and Nval > nlim
           if balise :
             # 0 SNPBim, 1 SNPSumStat, 2 CHR, 3 BP, 4 A1, 5 A2, 6 Z, 7 BETA, 8 SE, 9 AF, 10 N, P 11 
             rsbim="NA"
             if balisersbim==False:
                keybim='NA'
             else :
                 keybim=dickeyrs[key]
             dicres[key]=[keybim,gv(splline, rs_headerp, key), gv(splline, chr_headerp), gv(splline, bpval_headerp), gv(splline, a1_headerp).upper(), gv(splline, a2_headerp).upper(), gv(splline, z_headerp,formatfct=float), gv(splline, beta_headerp,formatfct=float), gv(splline, se_headerp,formatfct=float), gv(splline,af_headerp,formatfct=float), gv(splline, n_headerp, n_value_default,formatfct=float), gv(splline, pval_headerp, 'NA',formatfct=float)]  
             listnewkey.add(key)
    return (dicres, listnewkey)

def ExtractFreqN(bfile, dicsumstat, listkey,dicbim,freq_header,rs_header, n_header, chr_header,bpval_header, bin_plk, keep, threads, memory) :
   if (n_header or nvalue) and (freq_header is not None):
       return (dicsumstat, listkey)
   plkfreqfil=os.path.basename(bfile)
   out_range="tmp.frq"
   if bfile is None :
     print("no header for n or freq and bfile")
     sys.exit()
   Cmd=bin_plk+" -bfile "+bfile+" --freq --keep-allele-order "
   rangelist=out_range+".bed"
   writebed=open(rangelist, 'w') 
   [writebed.write("\t".join([dicsumstat[key][2],dicsumstat[key][3], dicsumstat[key][4], dicsumstat[key][0]])+'\n') for key in listkey]
   writebed.close()
   if keep is not None:
     Cmd+=" --keep "+keep
   Cmd+=" --threads "+str(threads)
   Cmd+=" --extract range "+rangelist
   Cmd+=" --memory "+memory
   Cmd+=" --out "+plkfreqfil
   os.system(Cmd)
   data_n=pd.read_csv(plkfreqfil+".frq",delim_whitespace=True)
   readfreq=open(plkfreqfil+".frq")
   #CHR	Chromosome code
   #SNP	Variant identifier
   #A1	Allele 1 (usually minor)
   #A2	Allele 2 (usually major)
  #MAF	Allele 1 frequency
  #NCHROBS	Number of allele observations
   header=readfreq.replace('\n', '').split()
   posrs=0
   posa1=1
   posa2=2
   posaf=3
   posn=4
   if maf is not None:
      mafupper= 1 - maf
   newlistkey=set([])
   for line in readfreq :
      splline=line.replace('\n', '').split()
      key=dicbim[splline[0]]
      # 0 SNP1, 1 SNP2, 2 CHR, 3 BP, 4 A1, 5 A2, 6 Z, 7 BETA, 8 SE, 9 AF, 10 N, P 11 
      if not freq_header : 
        af=float(splline[posaf])
        balise=True
        if maf is not None:
           balise=af>=maf and af<=mafupper
        if splline[posa1]==dicsumstat[key][5] and splline[posa2]==dicsumstat[key][4] :
          af = 1 - af
        if balise :
          dicsumstat[key][9]=af
          newlistkey.add(key)
        else :
          print("key "+key)
          del dicsumstat[key]
      if not n_header : 
         dicsumstat[key][10]=int(splline[posn])/2
   return (dicsumstat,newlistkey)

def addz(dicsumstat, listkey ,z_header) :
 # 0 SNP1, 1 SNP2, 2 CHR, 3 BP, 4 A1, 5 A2, 6 Z, 7 BETA, 8 SE, 9 AF, 10 N, P 11 
  if z_header :
     return dicsumstat
  for key in listkey :
      if dicsumstat[key][7]!='NA' and dicsumstat[key][8]!='NA':
         dicsumstat[key][6]=float(dicsumstat[key][7])/float(dicsumstat[key][8])
      else :
         dicsumstat[key][6]='NA'
  return dicsumstat
 
def formatval(x, fct) :
  if x =='NA':
    return x
  return fct(x)

def addbetase(dicsumstat, listkey ,beta_header, se_header) :
   if beta_header and se_header :
       return dicsumstat
   for key in listkey :
    z, p, n=formatval(dicsumstat[key][6], float),formatval(dicsumstat[key][9], float),formatval(dicsumstat[key][10], int)
    if z!='NA' and p!='NA' and n!='NA':
     b=z/sqrt(2*p*(1-p)*(n+z**2))
     se=1/sqrt(2*p*(1-p)*(n+z**2))
     dicsumstat[key][7]=b
     dicsumstat[key][8]=se
    else :
      print("pos "+key+" : z "+str(z)+" p : "+str(p)+" n "+str(n)+" is na")
   return(dicsumstat)

def compute_usingp(dicsumstat, listkey,used_p, beta_header, se_header, z_header) :
  if args.used_p==0 :
    return(dicsumstat)
  for key in listkey :
    if dicsumstat[key][7]!='NA':
      if float(dicsumstat[key][7])> 0 :
        sens=1
      else :
        sens = -1 
    elif dicsumstat[key][6]!='NA' :
      if float(dicsumstat[key][6]) > 0 :
        sens = 1
      else :
        sens = -1
    else :
      print("oos "+ key+" z and beta NA")
      continue
    dicsumstat[key][6]=abs(stats.norm.ppf(1-dicsumstat[key][11]/2))*sens
    z, p, n=float(dicsumstat[key][6]),float(dicsumstat[key][9]),int(dicsumstat[key][10])
    b=z/sqrt(2*p*(1-p)*(n+z**2))
    se=1/sqrt(2*p*(1-p)*(n+z**2))
    dicsumstat[key][7]=b
    dicsumstat[key][8]=se
  return dicsumstat

def checkrs(dicsumstat, listkey) :
   for key in listkey :
      #0 SNP1, 1 SNP2, 2 CHR, 3 BP, 4 A1, 5 A2, 6 Z, 7 BETA, 8 SE, 9 AF, 10 N, P 11 
      if dicsumstat[key][0]=='NA' :
         dicsumstat[key][0]=dicsumstat[key][1]
   return dicsumstat




#    plink_format.py --inp_asso $filegwas --chro_header ${params.head_chr_sumstat1} --bpval_header ${params.head_bp_sumstat1} --a1_header ${params.head_A1_sumstat1} --a2_header ${params.head_A2_sumstat1}   --pval_header ${params.head_pval_sumstat1}  --out ${out} --rs_header ${params.head_rs_sumstat1} $se $beta $z


def parseArguments():
    parser = argparse.ArgumentParser(description='extract rs in gwas file')
    parser.add_argument('--inp_asso',type=str,required=True)
    parser.add_argument('--rs',type=str,required=False, help="list rs", default="")
    parser.add_argument('--chro_header',type=str,required=True,help="chro header in inp files")
    parser.add_argument('--n_header',type=str,required=False,help="chro header in inp files")
    parser.add_argument('--bp_header',type=str,required=True,help="pos header in inp files")
    parser.add_argument('--beta_header',type=str,required=False,help="beta header in inp files")
    parser.add_argument('--z_header',type=str,required=False,help="Z header in inp files")
    parser.add_argument('--se_header',type=str,required=False,help="se header in inp files")
    parser.add_argument('--pval_header',type=str,required=True,help="p-value header in inp files")
    parser.add_argument('--rs_header',type=str,required=False,help="beta header in inp files")
    parser.add_argument('--n_lim',type=int,required=False,help="beta header in inp files", default=0)
    parser.add_argument('--a1_header',type=str,required=True,help="beta header in inp files")
    parser.add_argument('--keep_genet',type=str,help="beta header in inp files", default=0)
    parser.add_argument('--a2_header',type=str,required=True,help="beta header in inp files")
    parser.add_argument('--freq_header',type=str,required=False,help="frequencies header in inp files")
    parser.add_argument('--maf',type=float,default=0.0,help="minor allele frequencies")
    parser.add_argument('--bin_plk',type=str,required=False,help="plink binary", default="plink")
    parser.add_argument('--bfile',type=str,required=False,help="bfile if need to compute frequency or N", default=None)
    parser.add_argument('--keep',type=str,required=False,help="file of data used for if need to compute frequency or N", default=None)
    parser.add_argument('--used_p',type=int,required=False,help="file of data used for if need to compute frequency or N", default=0)
    parser.add_argument('--threads',type=int,required=False,help="", default=1)
    parser.add_argument('--memory',type=int,required=False,help="", default="1000")
    parser.add_argument('--n',required=False, help="bim file ")
    parser.add_argument('--out',type=str,default="out",help="b)")
    args = parser.parse_args()
    return args

args = parseArguments()

## keep information bim bam
print(" reading bfile :begin")
#(listkey, listdup, dicrskey,dickeyrs)
(listkey, listdup, dicrskey,dickeyrs)=extractinfobim(args.bfile+".bim")
print(len(listkey))
print(" reading bfile : end")
print(" reading sumstat: begin")
#def readsummarystat(sumstat, dickeyrs,listdup, listkey, rs_header,n_header, chr_header, bpval_header, beta_header, z_header, a1_header, a2_header, se_header, af_header, pval_header,n_value, maf) :
(dicsumstat,listkey)=readsummarystat(args.inp_asso, dickeyrs,listdup, listkey, args.rs_header, args.n_header, args.chro_header, args.bp_header, args.beta_header, args.z_header, args.a1_header, args.a2_header, args.se_header, args.freq_header, args.pval_header,args.n, args.maf, args.n_lim ,args.keep_genet) 
print(len(listkey))
print(" reading sumstat: end")
print(" reading add N and freq: begin")
(dicsumstat,listkey)=ExtractFreqN(args.bfile, dicsumstat, listkey,dicrskey,args.freq_header,args.rs_header, args.n_header, args.chro_header,args.bp_header, args.bin_plk, args.keep, args.threads,args.memory)
print(len(listkey))
print(" reading add N and freq: End")
print(" reading add beta se")
dicsumstat=addbetase(dicsumstat, listkey ,args.beta_header, args.se_header)
print(" reading add beta se :end")
dicsumstat=addz(dicsumstat, listkey ,args.z_header)
dicsumstat=compute_usingp(dicsumstat, listkey,args.used_p, args.beta_header, args.se_header, args.z_header)
dicsumstat=checkrs(dicsumstat, listkey)

writeplink=open(args.out+".plink",'w')
#writelz=open(args.out+".lz",'w')
#writegcta=open(args.out+".gcta",'w')
#0 SNP1, 1 SNP2, 2 CHR, 3 BP, 4 A1, 5 A2, 6 Z, 7 BETA, 8 SE, 9 AF, 10 N, P 11 
# 0 SNPBim, 1 SNPSumStat, 2 CHR, 3 BP, 4 A1, 5 A2, 6 Z, 7 BETA, 8 SE, 9 AF, 10 N, P 11 
headplink=["SNPKEY","SNP","SNPSS","CHR","BP","A1","A2", "Z","BETA","SE", "FRQ","N","P"]
headgcta=["SNPKEY","SNP","SNPSS","chr","bp","A1","A2", "z","b","se", "freq","N","p"]
#["#CHROM","BEGIN","END","MARKER_ID","PVALUE"]
headlz=["#CHROM","BEGIN","END","MARKER_ID","PVALUE"]

writeplink.write("\t".join(headplink)+"\n")
#writegcta.write("\t".join(headgcta)+"\n")

#writelz.write("\t".join(headlz)+"\n")
for key in listkey :
    info=dicsumstat[key]
    strdata=[str(x) for x in info]
    writeplink.write(key+"\t"+"\t".join(strdata)+"\n")
    #writegcta.write(key+"\t"+"\t".join(strdata)+"\n")
    #writelz.write("\t".join([strdata[2],strdata[3],strdata[3],strdata[1],strdata[11]])+"\n")




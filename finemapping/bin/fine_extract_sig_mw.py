#!/usr/bin/env python3

# Converts PLINK covariate and fam file into a covariate file for Gemma

import sys
import pandas as pd
import argparse
import numpy as np
import os
import scipy.stats as stats

def readbim(bimfile, sumstat) :
   readbim=open(bimfile)
   dicchro_rs={}
   dicchro_pos={}
   listsave=set([])
   #1       1:713250:G:C    0       713250  C       G
   for line in readbim :
      [chro,rs,cm,pos,A1,A2]=line.replace('\n','').split() 
      pos=int(pos)
      if (chro in sumstat) and (pos in sumstat[chro]):
       key=chro+':'+str(pos)
       if key in listsave :
         del dicchro_rs[chro][dicchro_pos[chro][pos][0]]
         del dicchro_pos[chro][pos]
         del sumstat[chro][pos]
         continue
       listsave.add(key)
       if (A1==  sumstat[chro][pos][0] and A2==sumstat[chro][pos][1]) or (A1==  sumstat[chro][pos][1] and A2==sumstat[chro][pos][0]) :
        sumstat[chro][pos].append(rs)
        if chro not in dicchro_rs :
           dicchro_rs[chro]={}
           dicchro_pos[chro]={}
        dicchro_rs[chro][rs]=pos
        dicchro_pos[chro][pos]=[rs,A1,A2]
       else :
        print(A1+"!="+sumstat[chro][pos][0])
   for chro in sumstat : 
      listbp=list(sumstat[chro].keys())
      for bp in listbp:
        if bp not in dicchro_pos[chro] :
          del sumstat[chro][bp]
   return (dicchro_pos,dicchro_rs,sumstat)
       
def extractinfoclump(clumpfile, wind) :
 readclump=open(clumpfile)
 dicrs={}
 headcl=readclump.readline()
 for line in readclump :
     spll=line.replace('\n','').split()
     if len(spll)<3 :
       continue
     chro=spll[0]
     rs=spll[2]
     bp=int(spll[3])
     listrs=spll[11]
     if listrs=='NONE' :
        listrs=[]
     else :
        listrs=listrs.replace('(1)','').replace('(2)','').split(',')
     listrs.append(rs)
     if chro not in dicrs :
         dicrs[chro]={}  
     dicrs[chro][bp]=[rs, listrs, bp-wind*1.5, bp+wind*1.5]
 return dicrs

def defined_wind(clumpres,bimfile_byrs, wind):
  listwind=[]
  for chro in clumpres.keys() :
     listbp=list(clumpres[chro].keys())
     for bp in listbp:
      infobp=clumpres[chro][bp]
      listpos=[]
      for x in infobp[1] :
        if x in bimfile_byrs[chro] :
           listpos.append(bimfile_byrs[chro][x]) 
      if len(listpos)==0 :
        print("warning "+str(chro)+" bp "+str(bp)+" deleted")
        del clumpres[chro][bp]
      else :
        posmin=min(listpos)
        posmax=max(listpos)
        clumpres[chro][bp].append(posmin)  
        clumpres[chro][bp].append(posmax)  
        clumpres[chro][bp][2]=min([posmin,bp-wind])
        clumpres[chro][bp][3]=max([posmax,bp+wind])
  return clumpres


def read_sumstat(filesumstat, clumpres, wind, chro_header, pos_header, a1_header, a2_header,beta_header, se_header, z_header,n_header ,freq_header, p_header, nval) :
  def getval(listval,posindex, Type=None):
     if posindex :
        val=listval[posindex] 
        if Type :
          if Type=='p' :
            if val.lower() == 'na':
               return None
            val = float(val)
            if val >1 or val<0 :
                  print(' val from ', p_header, ' >1 or < 0', str(val))
                  exit(2)
        return val    
     return None
  wind=wind*1.5
  readsumstat=open(filesumstat)
  header=readsumstat.readline().replace('\n','').split()
  chrohead=header.index(chro_header)
  bphead=header.index(pos_header)
  a1head=header.index(a1_header)
  a2head=header.index(a2_header)
  phead=header.index(p_header)
  nbcol=len(header)
  betahead=None
  if beta_header:
   betahead=header.index(beta_header)
  sehead=None
  if se_header:
   sehead=header.index(se_header)
  zhead=None
  if z_header:
   zhead=header.index(z_header)
  nhead=None
  if n_header:
   nhead=header.index(n_header)
  afhead=None
  if freq_header:
   afhead=header.index(freq_header)
  sumstat={}
  listchroclump=clumpres.keys()
  listsave=set([])
  for line in readsumstat :
    spl=line.replace('\n','').split()
    if len(spl)!=nbcol :
       print("line "+line)
       print("warning line doesn't have good format contains "+str(len(spl))+ " column and header "+str(nbcol))
       continue
    if spl[chrohead] not in listchroclump :
      continue
    key=spl[chrohead]+':'+spl[bphead]
    pos=int(spl[bphead])
    chro=spl[chrohead]
    if key in listsave :
      print("del "+chro+' '+str(pos))
      del sumstat[chro][pos]
      continue
    listsave.add(key)
    for posclmp in clumpres[chro].keys():
      if pos>=clumpres[chro][posclmp][2] and pos<=clumpres[chro][posclmp][3] :
        if chro not in sumstat :
          sumstat[chro]={}
        if nhead :
          n=getval(spl,nhead)
        elif nval:
          n = nval 
        else :
          n = None
        #          #['T', 'C', -0.4780229, 1.919196, -0.2490745603888295, 2010.0, 0.08209, 25000.0, '2:45922143:C:T']
        sumstat[chro][pos]=[getval(spl,a1head), getval(spl,a2head), getval(spl,betahead), getval(spl,sehead),getval(spl,zhead) , n, getval(spl,afhead, 'p'), getval(spl,phead,  'p')]
  return sumstat
#def extractinfobim(chro, begin, end,bimfile):
#   readbim=open(bimfile) 
#   listpos=[]
#   listref=[]
#   listalt=[]
#   listsnp=[]
#   listdup=set([])
#   for line in readbim :
#      splline=line.split()
#      pos=int(splline[3])
#      if chro == splline[0] and pos>=begin and pos<=end and (pos not in listdup):
#         if pos in listpos :
#           posindex=listpos.index(pos)
#           listpos.pop(posindex) 
#           listref.pop(posindex) 
#           listalt.pop(posindex) 
#           listdup.add(pos)
#           listsnp.pop(posindex)
#         else : 
#          listpos.append(pos)
#          listsnp.append(splline[1])
#          listref.append(splline[4])
#          listalt.append(splline[5])
#   return (listpos, listref, listalt, listsnp)

def appendfreq(bfile, result,biminfo, freq_header,rs_header, n_header, nval, bin_plk, keep, threads) :
   if freq_header and (n_header or nval):
        return result 
   out_range="range.bed"
   writebed=open(out_range, 'w')
   plkfreqfil=os.path.basename(bfile)
   for chro in result :
     for bp in result[chro] :
       writebed.write(chro+'\t'+str(bp)+'\t'+str(bp)+'\t'+chro+':'+str(bp)+'\n')
   writebed.close()
   if bfile==None :
     print("no header for n or freq and bfile")
     sys.exit()
   Cmd=bin_plk+" -bfile "+bfile+" --freq --keep-allele-order --extract  range "+out_range
   if keep :
     Cmd+=" --keep "+keep
   Cmd+=" --threads "+str(threads)
   Cmd+=" --out "+plkfreqfil
   if os.path.exists(plkfreqfil+".frq") ==False :
     print(Cmd)
     os.system(Cmd)
   readfreq=open(plkfreqfil+".frq")
   splheader=readfreq.readline().replace('\n','').split()
   poschro=splheader.index("CHR") 
   posrs=splheader.index("SNP") 
   posN=splheader.index("NCHROBS")
   posA1=splheader.index("A1")
   posA2=splheader.index("A2")
   posMAF=splheader.index("MAF")
   #[getval(spl,a1head), getval(spl,a2head), getval(spl,betahead), getval(spl,sehead),getval(spl,zhead) , getval(spl,nhead), getval(spl,afhead)]
   for line in readfreq :
       spll=line.replace('\n','').split()
       chro=spll[poschro]
       rs=spll[posrs]
       try :
         bp=biminfo[chro][rs]
       except :
         print("rs"+rs+" not found in bim")
         exit(-1)
       if n_header==None :
          result[chro][bp][5]=round(int(spll[posN])/2)
       if freq_header==None :
         try :
            maf=float(spll[posMAF]) 
            A1=spll[posA1]
            A2=spll[posA2]
            if A1 != result[chro][bp][0] :
              maf=1 - maf 
            result[chro][bp][6] = maf
         except :
            print("maf error "+ line+' '+spll[posMAF])
   for chro in result :
     listbp=list(result[chro].keys())
     for bp in listbp :
       if (not result[chro][bp][5])  or (not result[chro][bp][6]):
          print("warning "+str(chro)+" bp "+str(bp)+" deleted")
          del result[chro][bp]
   return result

def checksumstat(sumstat,maf) :
   listchro=sumstat.keys()
   for chro in listchro: 
      listbp_chro=list(sumstat[chro].keys())
      for bp in listbp_chro:
         #[getval(spl,a1head), getval(spl,a2head), getval(spl,betahead), getval(spl,sehead),getval(spl,zhead) , getval(spl,nhead), getval(spl,afhead)]
         [a1,a2,beta, se, z, n, af, p, rsbim]=sumstat[chro][bp]
         if not p :
           del sumstat[chro][bp]
           continue 
         n = float(n)
         p = float(p)
         if args.z_pval==1 :
          z=stats.norm.ppf(1-p/2)
          if z < 0 :
             z = - z
          if float(beta) < 0:
             z = - z
             beta = z * float(se)
         elif z==None :
          z = float(beta)/ float(se)
         if beta == None :
             beta=z/sqrt(2*p*(1-p)*(n+z**2))
             se=1/sqrt(2*p*(1-p)*(n+z**2))
         beta=float(beta)
         se=float(se)
         af=float(af) 
         if af > 0.5 :
           af= 1 - af
           beta =  - beta
           z =  - z
           a1tmp=a1
           a1=a2
           a2=a1tmp
         if af>=maf :
           sumstat[chro][bp]=[a1,a2,beta, se, z, n, af, p,rsbim]
         else :
           print("maf del :"+chro+' '+str(bp))
           del sumstat[chro][bp]
   return sumstat 
EOL=chr(10)

def parseArguments():
    parser = argparse.ArgumentParser(description='extract rs in gwas file')
    parser.add_argument('--inp_resgwas',type=str,required=True)
    parser.add_argument('--clump',type=str,required=True)
    parser.add_argument('--min_pval',type=float,required=False, help="list rs")
    parser.add_argument('--chro_header',type=str,required=True,help="chro header in inp files")
    parser.add_argument('--n_header',type=str,required=False,help="chro header in inp files")
    parser.add_argument('--pos_header',type=str,required=True,help="pos header in inp files")
    parser.add_argument('--beta_header',type=str,required=True,help="beta header in inp files")
    parser.add_argument('--se_header',type=str,required=True,help="beta header in inp files")
    parser.add_argument('--p_header',type=str,required=True,help="beta header in inp files")
    parser.add_argument('--rs_header',type=str,required=False,help="beta header in inp files")
    parser.add_argument('--a1_header',type=str,required=True,help="beta header in inp files")
    parser.add_argument('--z_header',type=str,required=False,help="beta header in inp files")
    parser.add_argument('--around',type=float,required=False,help="beta header in inp files", default=50000)
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
clumpres=extractinfoclump(args.clump, args.around)
print('clump',' '.join([x+':'+str(len(clumpres[x])) for x in clumpres]))
sumstat=read_sumstat(args.inp_resgwas, clumpres,  args.around, args.chro_header, args.pos_header, args.a1_header, args.a2_header, args.beta_header, args.se_header, args.z_header,args.n_header ,args.freq_header, args.p_header, args.n)
print('sumstat',' '.join([x+':'+str(len(sumstat[x])) for x in sumstat]))
(diccbim_pos,diccbim_rs, sumstat)=readbim(args.bfile+'.bim', sumstat)
print('afterbim',' '.join([x+':'+str(len(sumstat[x])) for x in sumstat]))
sumstat=appendfreq(args.bfile, sumstat,diccbim_rs, args.freq_header,args.rs_header, args.n_header, args.n,args.bin_plk, args.keep, args.threads)
clumpres=defined_wind(clumpres,diccbim_rs, args.around)
sumstat=checksumstat(sumstat, args.maf)

## writing      
#0 1 2 3 4 
#[a1 0 ,a2 1 ,beta 2 , se 3 , z 4 , n 5 , af 6 , rsbim 7 ,p 8]

#posgcta=[7,0,1,6,2, 3,8, 5]
headgcta=['SNP','A1','A2','freq','b','se','p','N']
headfm=["rsid","chromosome","position","allele1","allele2","maf", "beta", "se"]
#"rsid","chromosome","position","allele1","allele2","maf","beta","se", "p","z"
headall=["rsid","chromosome","position","allele1","allele2","maf", "beta", "se", "z", 'p', 'n']
##out_fileZ=args.out_head+"_caviar.z"
#out_all=args.out_head+".all"
#small.to_csv(out_all, sep=TAB, header=True, index=False,na_rep="NA")
posfinemap=[7,]
maxn=0
for chro in clumpres :
   for bp in clumpres[chro] :
      key=chro+'_'+str(bp)
      out_filefm=key+"_finemap.z"
      out_filecaviar=key+"_caviar.z"
      out_range=key+".range"
      out_gcta=key+'.gcta'
      out_all=key+'.all'
      out_paintor=key+'.paintor'
      out_pos=key+'.pos'
      out_rs=key+'.rs'
      begin=clumpres[chro][bp][2]
      end=clumpres[chro][bp][3]
      listsumstat=sumstat[chro].keys() 
      listsumstat=[x for x in listsumstat if x >= begin and x<= end]
      listsumstat.sort()
      if len(listsumstat)>0 :
        writegcta=open(out_gcta ,'w')
        writefm=open(out_filefm,'w')
        writecaviar=open(out_filecaviar,'w')
        writerange=open(out_range,'w')
        writeall=open(out_all, 'w')
        writepaintor=open(out_paintor,'w')
        writepos=open(out_pos,'w')
        writers=open(out_rs,'w')
        writegcta.write(' '.join(headgcta)+'\n')
        writefm.write(' '.join(headfm)+'\n')
        writeall.write(' '.join(headall)+'\n')
        writepaintor.write('Z\n')
        writepos.write("chromosome position\n")
        for bpsumstat in listsumstat :
          [a1,a2,beta , se , z , n , af , p, rsbim ]=sumstat[chro][bpsumstat]
          maxn=max(n,maxn)
          writegcta.write(' '.join([str(x) for x in [rsbim, a1,a2,af,beta,se,p,n]])+'\n')
          writefm.write(' '.join([str(x) for x in [rsbim,chro,bpsumstat, a1,a2,af,beta,se]])+'\n')
          writecaviar.write(" ".join([str(x) for x in [rsbim,z]])+'\n')
          writerange.write(' '.join([str(x) for x in [chro,bpsumstat, bpsumstat,rsbim]])+'\n')
          writeall.write(' '.join([str(x) for x in [rsbim,chro,bpsumstat,a1,a2,af,beta,se,z,p,n ]])+'\n')
          writepaintor.write(str(z)+'\n')
          writepos.write(chro+' '+str(bpsumstat)+'\n')
          writers.write(rsbim+'\n')
        writegcta.close()
        writefm.close()
        writecaviar.close()
        writerange.close()
        writeall.close()
        writepaintor.close()
        writepos.close()
        writers.close()
      
writen=open('n.out', 'w')
writen.write(str(maxn))
writen.close()
        


#result[args.chro_header]=result[args.chro_header].astype(str)
#result.drop_duplicates(subset=[args.chro_header, args.pos_header],keep=False,inplace=True)
#rs=args.
#freq_header=args.freq_header
#n_header=args.n_header
#rs_header=args.rs_header
#maf=args.maf
#out_file=args.out_head+"_finemap.z"
#out_fileZ=args.out_head+"_caviar.z"
#if args.rs :
#   sub_result=result[result[args.rs_header]==rs]
#   chro=sub_result[args.chro_header].tolist()[0]
#   pos=sub_result[args.pos_header].tolist()[0]
#   begin=pos-args.around
#   end=pos+args.around
#else :
#   chro=args.chro
#   begin=args.begin
#   end=args.end
#TAB=chr(9)
###rsid chromosome position allele1 allele2 maf beta se
#
#for rs in listrs :
#  (listbim, listalt, listref, listsnp)=extractinfobim(chro, begin,end,args.bfile+".bim")
#  small=result[(result[args.chro_header]==chro) & (result[args.pos_header]>=begin) & (result[args.pos_header]<=end)]
#  small=pd.merge(small, pd.DataFrame({"pos":listbim, 'altbim232':listalt, 'refbim232':listref, 'snpbim232':listsnp}),left_on=args.pos_header, right_on="pos" )
#small=small[((small['altbim232']==small[args.a1_header]) & (small['refbim232']==small[args.a2_header])) | ((small['altbim232']==small[args.a2_header]) & (small['refbim232']==small[args.a1_header]))]
#
#if args.min_pval and small[args.p_header].min()> args.min_pval:
#   sys.exit('min pval '+str(args.min_pval)+'>'+' min(p) '+ str(small[args.p_header].min()))
#small=small.loc[small[args.pos_header].isin(listbim)]
#if args.n : 
#  small['N']=args.n
#  n_header='N'
#
#if (n_header==None or freq_header==None):
#  (freq_header, n_header, small)=appendfreq(args.bfile, small, freq_header,rs_header, n_header, args.chro_header,args.pos_header,args.bin_plk, args.keep, args.threads)
#print(small)
#
#PosCol=['snpbim232',args.chro_header, args.pos_header, args.a1_header, args.a2_header, freq_header, args.beta_header, args.se_header, args.p_header, n_header]
#RsName=["rsid","chromosome","position","allele1","allele2","maf", "beta", "se", "p","N"]
#small=small[(small[rs_header]==rs) | ((small[freq_header]>maf) & (small[freq_header]<(1-maf)))]
#
#small=small[PosCol] 
#small.columns=RsName
#small=small.sort_values(["position"])
#### for gcta
#smallgcta=small[["rsid","chromosome","position","allele1","allele2","maf", "beta", "se", 'p','N']]
#smallgcta=smallgcta.rename(columns={"rsid": "SNP", "chromosome": "chr", "position":'bp', 'allele1':'A1', 'allele2':'A2', 'beta':'b', 'maf':'freq'})
#
#out_gcta=args.out_head+'.gcta'
#
### change freq and allele
#bal=small['maf']>0.5
#small['allele1_tmp']=small['allele1']
#small['allele2_tmp']=small['allele2']
#small['maf_tmp']=small['maf']
#small['beta_tmp']=small['beta']
#small.loc[bal,'allele1']=small.loc[bal,'allele2_tmp']
#small.loc[bal,'allele2']=small.loc[bal,'allele1_tmp']
#small.loc[bal,'beta']= - small.loc[bal,'beta_tmp']
#small.loc[bal,'maf']= 1 - small.loc[bal,'maf_tmp']
#
#if args.z_pval==0 :
#  small['Z']=small['beta']/small['se']
#else :
#  tmpbeta=small['beta'].copy()
#  tmpbeta[tmpbeta>=0]=1
#  tmpbeta[tmpbeta<=0]=-1
#  small['Z']=stats.norm.ppf(1-small['p']/2)
#  small['Z']=small['Z'].abs()*tmpbeta
#  small['beta']=small['Z']*small['se']
#
#
#
#out_range=args.out_head+".range"
#
#out_all=args.out_head+".all"
#
#out_range=args.out_head+".paintor"
#out_pos=args.out_head+'.pos'
#small[["chromosome","position"]].to_csv(out_pos, sep=" ", header=True, index=False,na_rep="NA")
#
#out_pos=args.out_head+'.rs'
#small[["rsid"]].to_csv(out_pos, sep=" ", header=True, index=False,na_rep="NA")
#
#

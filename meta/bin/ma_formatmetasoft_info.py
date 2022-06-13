#!/usr/bin/env python3

########################################################
# ma_formatmetasoft.py modification of plink2metasoft.py                                    
#   Convert file of assoc files to Metasoft input file  
#   Free license -- you are free to use it in any ways 
#   Buhm Han (2012)                                    
#   memories used is very high
########################################################

import sys, subprocess, os
#import stats
from scipy import stats

comple = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
vectortorbase=['A','T','C','G']

rsHead="RSID"
ChHead="CHRO"
PsHead="POS"
A1Head="A1"
A2Head="A2"
BetHead="BETA"
PHead="PVAL"
SError="SE"
NHead="N"
FreqHead="FREQA1"

# PROCESS ARGUMENTS
if len(sys.argv) < 3:
    print_and_exit()
out=sys.argv[1]
files=sys.argv[2:]

# MERGE STUDIES
fout=open(out+'.meta','w')
foutN=open(out+'.N','w')
#fmap=open(out+'.mmap','w')
flog=open(out+'.log','w')


def GetInfoRsGWAS(rsid, snp,A1Pivot, A2Pivot,CompSE, PosA1Head, PosA2Head,  PosBetHead,PosN, PosFreq) :
    snp[PosA1Head]=snp[PosA1Head].upper()
    snp[PosA2Head]=snp[PosA2Head].upper()
    studyindex=-1
    pivotstudyindex=-2
    beta=snp[PosBetHead]
    if beta == "NA" :
        return (False,None, None)
    try :
      float(beta)
    except :
      return (False,None, None)
    if PosA1Head : 
      if snp[PosA1Head] not in vectortorbase:
         return (True,None, None)
    if PosA2Head : 
      if snp[PosA2Head] not in vectortorbase:
         return (True,None, None)
    if PosN :  
       NVal=int(snp[PosN])
    else :
       NVal=None
    if PosFreq :  
       try :
        FreqVal=float(snp[PosFreq])
       except :
        FreqVal=None
    else :
       FreqVal=None
    # CHECK ALLELE TO PIVOT
    if PosA2Head and PosA1Head :
       if A1Pivot == snp[PosA1Head] and A2Pivot == snp[PosA2Head]:
          return (True,NVal, FreqVal)
       elif A1Pivot == snp[PosA2Head] and A2Pivot == snp[PosA1Head]:
            # SIMPLE FLIP
            FreqVal = 1 - FreqVal
       elif A1Pivot == comple[snp[PosA1Head]] and A2Pivot == comple[snp[PosA2Head]]:
                        # STRAND INCONSIS., BUT GOOD
            flog.write('FLIP_STRAND %s in study %d\n'%(rsid,studyindex))
       elif  A1Pivot == comple[snp[PosA2Head]] and A2Pivot == comple[snp[PosA1Head]]:
             # STRAND INCONSIS., SIMPLE FLIP
             flog.write('FLIP_STRAND %s in study %d\n'%(rsid,studyindex))
             if FreqVal :
                FreqVal = 1 - FreqVal
             else :
               return (True,None, None)
       else:
            return (True,None, None)
    return (True,NVal, FreqVal)


# READ FILES
#studies=[]
newfileslist=[]
rsidschar={}
rsidsN={}
rsidsFreq={}
rsidsinfo={}
CmtFile=0
listrsall=set([])

for f in files:
    listrsfile=set([])
    study={}
    fin=open(f)
    colnames=fin.readline().replace('\n','').split()
    CompSE=False
    if SError not in colnames:
       CompSE=True
       if PHead not in colnames :
          print("no Phead in "+ f+" : skip")
          continue
    if BetHead not in colnames :
       print("no beta in "+ f+" : skip")
       continue
    newfileslist.append(os.path.basename(f))
    PosA1Head=-1
    PosA2Head=-1
    PosPHead=None
    PosSError=None
    PosPosHead=-1
    PosChroHead=-1
    if A1Head in colnames :
       PosA1Head=colnames.index(A1Head)
    if A2Head in colnames:
       PosA2Head=colnames.index(A2Head)
    if PHead in  colnames:
       PosPHead=colnames.index(PHead)
    if SError in  colnames:
       PosSError=colnames.index(SError)
    if PsHead in colnames :
       PosPosHead=colnames.index(PsHead)
    if ChHead in colnames :
       PosChroHead=colnames.index(ChHead)
    if NHead in colnames :
       PosNHead=colnames.index(NHead)
    else :
       PosNHead=None
    if FreqHead in colnames :
       PosFreq=colnames.index(FreqHead)
    else :
       PosFreq=None
    PosBetHead=colnames.index(BetHead)
    PosRsHead=colnames.index(rsHead)
    for line in fin:
        spll=line.split()
        rsid=spll[PosRsHead]
        if rsid not in listrsall :
           listrsall.add(rsid)
           rsidsN[rsid]=[None]*CmtFile
           rsidsFreq[rsid]=[None]*CmtFile
           if PosA2Head>=0 and PosA1Head>=0:
             rsidsinfo[rsid]=[0,spll[PosA1Head],spll[PosA2Head]]

        if rsid not in listrsfile :
           if rsid in rsidsinfo :
             #def GetInfoRsGWAS(rsid, snp,A1Pivot, A2Pivot,CompSE, PosA1Head, PosA2Head,  PosBetHead,PosN, PosFreq) :
             (balise, N, Freq)=GetInfoRsGWAS(rsid, spll,rsidsinfo[rsid][1], rsidsinfo[rsid][2],CompSE, PosA1Head, PosA2Head, PosBetHead, PosNHead, PosFreq)
             if balise :
               rsidsinfo[rsid][0]+=1
             rsidsN[rsid].append(N)
             rsidsFreq[rsid].append(Freq)
             listrsfile.add(rsid)
        else :
           print("rs "+ rsid +" multi times :skip "+f)
    toappend=[x for x in listrsall if x not in listrsfile]
    for x in toappend :
      rsidsN[x].append(None)
      rsidsFreq[x].append(None)
    fin.close()
    CmtFile+=1

for rsid in rsidsN :
  if len(rsidsN[rsid])!=len(files):
     print(rsidsN[rsid])
     print(rsidsN[rsid])
     print("error len of "+rsid+" error") 
     exit(2)

def getsumandfreq(N, freq):
 cmt=0
 SumN=0
 SumN2=0.0
 SumAll=0.0
 StrN='' 
 StrFrq='' 
 for x in N :
    if x :
      SumN+=x
      if freq[cmt] :
        SumAll+=(freq[cmt]*x)
        SumN2+=x
        StrFrq+=str(freq[cmt])+"\t"    
      else  :
       StrFrq+="NA\t"
      StrN+=str(x)+"\t"    
    else :
     StrN+="NA\t"
     StrFrq+="NA\t"
    cmt+=1
 if SumN == 0 :
   NFinal=StrN+'NA'
   FreqFinal=StrFrq+'NA'
 elif SumN2 == 0 :
   NFinal=StrN+str(SumN)
   FreqFinal=StrFrq+'NA'
 else :
    NFinal=StrN+str(SumN)
    FreqFinal= StrFrq+str((SumAll/float(SumN2)))
 return (NFinal, FreqFinal)
     

for rsid in listrsall :
   if rsidsinfo[rsid][0] > 0:
         (NF,Freq)=getsumandfreq(rsidsN[rsid], rsidsFreq[rsid])
         foutN.write(rsid +"\t"+NF+"\t"+Freq+'\n')

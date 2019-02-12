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

# PROCESS ARGUMENTS
if len(sys.argv) < 3:
    print_and_exit()
out=sys.argv[1]
files=sys.argv[2:]

# READ FILES
studies=[]
newfileslist=[]
for f in files:
    study={}
    fin=open(f)
    colnames=fin.readline().split()
    if BetHead not in colnames :
       print("no beta in "+ f+" : skip")
       continue
    if PHead not in colnames :
       print("no Phead in "+ f+" : skip")
       continue
    newfileslist.append(os.path.basename(f))
    for line in fin:
        snp={}
        for (x,y) in zip(colnames,line.split()):
            snp[x]=y
        rsid=snp[rsHead]
        study[rsid]=snp
    fin.close()
    studies.append(study)
    
# UNION OF SNPS
allsnps={}
for study in studies:
    for rsid, snp in study.items():
        allsnps[rsid]=snp

# SORT SNPS
#rsids=sorted(allsnps.keys(), 
#             key=lambda x:(allsnps[x][ChHead],int(allsnps[x][PsHead])))

wf=open(out+'.'+'files', 'w') 
wf.write("\n".join(newfileslist))
wf.close()
rsids=allsnps

# MERGE STUDIES
fout=open(out+'.meta','w') 
#fmap=open(out+'.mmap','w')
flog=open(out+'.log','w')  
for rsid in rsids:
    output=rsid+'\t'
    pivot=-1
    pivotstudyindex=-1
    numstudy=0
    for study in studies:
        studyindex=studies.index(study)+1
        if rsid in study:
            snp=study[rsid]
            snp[A1Head]=snp[A1Head].upper()
            snp[A2Head]=snp[A2Head].upper()
            if snp[A1Head] not in vectortorbase  or snp[A1Head] not in vectortorbase :
               output+='NA NA '
               continue 
            if 'OR' in snp:
                beta='log(%s)'%snp['OR']
            elif BetHead in snp:
                beta=snp[BetHead]
            else:
                assert 0, 'OR or BETA must be in columns'
            if PHead in snp:
                p=snp[PHead]
            else:
                assert 0, 'P must be in columns'
            #if SError in snp :
            if SError in snp :
               stderr=snp[SError] 
            else :
               if float(p)==0.5 :
                  p2=0.499999
               else :
                  p2=p
               stderr="%f"%abs(float(beta)/stats.norm.ppf(float(p2), loc =0, scale = 1))  #'abs(%s/qnorm(%s/2))'%(beta,p)
            if pivot == -1:
                pivot=snp # 1ST STUDY's SNP INFO IS PIVOT
                pivotstudyindex=studyindex
            else:
                # CHECK ALLELE TO PIVOT
                if A2Head in pivot and A2Head in snp:
                    if pivot[A1Head] == snp[A1Head] and \
                       pivot[A2Head] == snp[A2Head]:
                        # GOOD
                        pass
                    elif pivot[A1Head] == snp[A2Head] and \
                         pivot[A2Head] == snp[A1Head]:
                        # SIMPLE FLIP
                        beta='%f'%(float(beta)*-1)
                    elif pivot[A1Head] == comple[snp[A1Head]] and \
                         pivot[A2Head] == comple[snp[A2Head]]:
                        # STRAND INCONSIS., BUT GOOD
                        flog.write('FLIP_STRAND %s in study %d\n'%(rsid,studyindex))
                    elif pivot[A1Head] == comple[snp[A2Head]] and \
                         pivot[A2Head] == comple[snp[A1Head]]:
                        # STRAND INCONSIS., SIMPLE FLIP
                        flog.write('FLIP_STRAND %s in study %d\n'%(rsid,studyindex))
                        beta='%f'%(float(beta)*-1)
                    else:
                        flog.write('EXCLUDE %s due to allele inconsistency: A1:%s A2:%s in study %d but A1:%s A2:%s in study %d\n'
                                   %(rsid, pivot[A1Head], pivot[A2Head], pivotstudyindex,
                                     snp[A1Head], snp[A2Head], studyindex))
                else: 
                    if pivot[A1Head] == snp[A1Head]:
                        # GOOD
                        pass
                    else:
                        flog.write('EXCLUDE %s due to allele inconsistency: A1:%s in study %d but A1:%s in study %d\n'
                                   %(rsid, pivot[A1Head], pivotstudyindex,
                                     snp[A1Head], studyindex))
                # CHECK CHR & BP TO PIVOT
                #if pivot[ChHead] != snp[ChHead]:
                #    flog.write('WARNING %s has chr inconsistency\n'%rsid)
                #if pivot[PsHead] != snp[PsHead]:
                #    flog.write('WARNING %s has basepair inconsistency\n'%rsid)
            output+=beta+' '+stderr+' '
            numstudy+=1
        else:
            output+='NA NA '
    if numstudy == 1:
        flog.write('EXCLUDE %s due to being in single study\n'%rsid)
    else:
        fout.write(output+'\n')
    #    fmap.write('%s\t%s\t%s\t%s\t%d\n'%
    #               (rsid, pivot[ChHead], pivot[PsHead], pivot[A1Head], numstudy))
fout.close()
#fmap.close()
flog.close()

# CALL R TO EVALUATE MATH EXPRESSION
#subprocess.call(['R --vanilla --slave "--args '+out+'.meta.tmp '+out+'.meta" < plink2metasoft_subroutine.R'],shell=True)
#subprocess.call(['rm '+out+'.meta.tmp'], shell=True)



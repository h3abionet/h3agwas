#!/usr/bin/env python3
import sys
import argparse

''' 
clean a vcf file and sumtat using a1 / a2
'''

def codeallele(chro,bp,a1,a2):
  if not a1 :
     return chro+':'+bp 
  if a1>a2 :
    return chro+':'+bp+':'+a1+':'+a2
  return chro+':'+bp+':'+a2+':'+a1

def read_annovar (fileannovar, baliseallele):
 def codeallele_anov(spl, baliseallele):
   if baliseallele :
      return (codeallele(spl[0],spl[2],spl[3],spl[4]),spl)
   else :
      return (codeallele(spl[0],spl[2],None,None),spl)
 read=open(fileannovar)
 tmp=[codeallele_anov(x.replace('\n','').split(), baliseallele) for x in read]
 read.close()
 dicout={}
 for (x,spl) in tmp :
  dicout[x]=spl
 return dicout



def read_bim(filebim, baliseallele)
 def codeallele_bim(spl, baliseallele):
   if baliseallele :
      return (codeallele(spl[0],spl[3],spl[4],spl[5]),spl)
   else :
      return (codeallele(spl[0],spl[2],None,None),spl)
 read=open(filebim)
 tmp=[codeallele_bim(x.replace('\n','').split(), baliseallele) for x in read]
 read.close()
 dicout={}
 for (x,spl) in tmp :
  dicout[x]=spl
 return dicout

def writebim (spl) :

def parseArguments():
    parser = argparse.ArgumentParser(description='merge file statistics with file frequence generate by plink')
    parser.add_argument('--gwas',type=str,required=True, help="association files")
    parser.add_argument('--vcf',type=str,help="vcf files")
    parser.add_argument('--out', type=str,help="",default="out")
    parser.add_argument('--baliseallele', type=int,help="",default=1)
    parser.add_argument('--format', type=str,help="",default='annovar')
    args = parser.parse_args()
    return args


args = parseArguments()
baliseallele=args.baliseallele
if args.format == 'annovar' :
  Format=args.format
  listkeyannov=read_annovar(args.gwas, baliseallele)
else args.format == 'bim' :
  listkeyannov=read_bim(args.gwas, baliseallele)
  baliseallele=True

keyannov=listkeyannov.keys()

# vcf contained new rs
readvcf = open(args.vcf)

if Format == 'bim' :
 writeabim=open(args.out+'.bim','w')
else  :
 writeabim=open(args.out+'.annov','w')

write_log=open(args.out+'.log','w')

for vcfl in readvcf :
  if vcfl[0] != "#" :
    spl=vcfl.replace('\n','').split()  
    if baliseallele :
      a2list=spl[4].split(',')
      for a2 in a2list:
        codeall=codeallele(spl[0],spl[1],spl[3],a2)
        if codeall in keyannov:
          writeabim.write("\t".join([spl[0],spl[1],spl[2],spl[3],a2, spl[2]])+'\n')
          spl[4]=a2
          write_vcf.write("\t".join(spl)+'\n')
          del listkeyannov[codeall]
    else :
       codeall=codeallele(spl[0],spl[1],None,None)
       if codeall not in keyannov:
          write_log.write('error '+codeall+' not found in '+args.gwas+'skip\n')
       else :
          write_vcf.write(vcfl)
          writeabim.write("\t".join([spl[0],spl[1],spl[1],spl[3],spl[4], spl[2]])+'\n')
          del listkeyannov[codeall]

for allele in listkeyannov :
   write_log.write('warning '+allele+' not found in '+args.vcf+'\n')
   if baliseallele :
       spl=listkeyannov[allele]
       writeabim.write('\t'.join([spl[0],spl[1],spl[2],spl[3],spl[4], spl[0]+':'+spl[1]]))

writeabim.close()
write_log.close()
    
    

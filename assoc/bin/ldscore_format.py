#!/usr/bin/env python3

import sys
import argparse
import gzip


""" 
transform ld score file to ensure that rsid is same with plink file
"""


EOL=chr(10)
def parseArguments():
    parser = argparse.ArgumentParser(description='transform ld score file to ensure that rsid is same with plink file')
    parser.add_argument('--ldscore',type=str,required=True)
    parser.add_argument('--plk_bim',type=str,required=True)
    parser.add_argument('--out',type=str,required=True)
    args = parser.parse_args()
    return args

args=parseArguments()
## extraction of bim file plink format, rsid by chr:pos
readplk=open(args.plk_bim)
dicplk={}
for lineplk in readplk :
    spllineplk=lineplk.split()
    if spllineplk[0] not in dicplk :
      dicplk[spllineplk[0]]={} 
    dicplk[spllineplk[0]][spllineplk[3]]=spllineplk[1]
readplk.close()

readldscore=gzip.open(args.ldscore, 'rb')
writeldscore=gzip.open(args.out, 'wb')
header=readldscore.readline()
writeldscore.write(header)
for lineld in  readldscore :
   splt=lineld.decode('utf-8').replace('\n','').split()
   #rs367896724     1       10177   2.435
   #keypos=splt[1]+':'+splt[2]
   chro=splt[1]
   bp=splt[2]
   if chro in dicplk and bp in dicplk[chro]:
      #print('replace '+splt[0]+' '+dicplk[chro][bp])
      splt[0]=dicplk[chro][bp]
      newline='\t'.join(splt)+'\n'
      writeldscore.write(newline.encode('utf-8'))
writeldscore.close()


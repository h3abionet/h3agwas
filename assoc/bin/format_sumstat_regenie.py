#!/usr/bin/env python3

# Converts PLINK covariate and fam file into a covariate file for Gemma
import sys
import argparse
import math

filesumstat=sys.argv[1]

readsum=open(filesumstat)
writesum=open(sys.argv[2], 'w')

headlin=readsum.readline().replace('\n', '')
splhead=headlin.split()
poslog10=splhead.index("LOG10P")
headlin=headlin+" P"+"\n"

writesum.write(headlin)

for line in readsum :
 line=line.replace('\n','')
 spli=line.split()
 try :
  P=str(math.exp(-float(spli[poslog10])))
 except :
  P="NA"
 writesum.write(line + " "+P+"\n")


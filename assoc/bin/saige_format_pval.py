#!/usr/bin/env python3
''' 
format file for gcta, append freq and n if need and bfile
'''
from  subprocess import CalledProcessError
import subprocess
import os
import argparse
import numpy as np
import pandas as pd
import sys

EOL = chr(10)
#SNP A1 A2 freq b se p N 
def parseArguments():
    parser = argparse.ArgumentParser(description='format file for gcta, append N and frequencie if not present using bed file')
    parser.add_argument('--inp_asso',type=str,required=True, help="association files")
    parser.add_argument('--out', type=str,help="out of format file",default="test.tex")
    args = parser.parse_args()
    return args

args = parseArguments()
inp     = args.inp_asso
out = args.out

readf=open(inp)
header=readf.readline().replace('\n','').split()
#  #CHR     POS     MarkerID        Allele1 Allele2 AC_Allele2      AF_Allele2      imputationInfo  BETA    SE      Tstat   var     p.value p.value.NA      Is.SPA  AF_case AF_ctrl N_case  N_ctrl
afcase=header.index('AF_case')
afcontrol=header.index('AF_ctrl')

ncase=header.index('N_case')
ncontrol=header.index('N_ctrl')

writegwas=open(out,'w')

header.append('N_total')
#header.append('af_total')
writegwas.write("\t".join(header)+'\n')
for line in readf :
  spl=line.replace('\n','').split()
  nbcase=float(spl[ncase])
  nbcontrol=float(spl[ncontrol])
  spl.append(str(int(nbcase+nbcontrol))) 
  #afcase=float(spl[afcase])
  #afcontrol=float(spl[afcontrol])
  #afall=(afcase*ncase + afcontrol*ncontrol)/(ncase + ncontrol)
  #spl.append(str(afall))
  writegwas.write("\t".join(spl)+'\n')

writegwas.close()
readf.close()


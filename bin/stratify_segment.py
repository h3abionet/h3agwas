#!/usr/bin/env python3

# Converts PLINK covariate and fam file into a covariate file for Gemma

import sys
import pandas as pd
import argparse
import numpy as np



EOL=chr(10)
def parseArguments():
    parser = argparse.ArgumentParser(description='fill in missing bim values')
    parser.add_argument('--inp_ld',type=str,required=True)
    parser.add_argument('--nb_split',type=float,required=True)
    parser.add_argument('--out',type=str,required=True)
    args = parser.parse_args()
    return args

args=parseArguments()
data = pd.read_csv(args.inp_ld,delim_whitespace=True)
splitlist=np.arange(0,1,1.0/float(args.nb_split)).tolist()
splitlist+=[1]
quant=np.quantile(data[["ldscore_region"]],splitlist)
quant[0]=quant[0]-1
quant[len(quant)-1]=quant[len(quant)-1]+1

for x in range(0,len(quant)-1) :
   listsnp=data[(data["ldscore_region"] >quant[x]) & (data["ldscore_region"] <=quant[x+1])]
   listsnp=listsnp["SNP"].tolist()
   ecrire=open(args.out+"_"+str(x)+".reg", "w") 
   ecrire.write("\n".join(listsnp))
   ecrire.close()



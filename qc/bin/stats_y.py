#!/usr/bin/env python3
#Load HWE P-value file and generate frequency_distribution
 
import pandas as   pd
import numpy  as np
import math
import sys
from matplotlib import use
use('Agg')
import argparse
import matplotlib.pyplot as plt
import sys
import scipy.stats as stats 
 
def parseArguments():
    parser = argparse.ArgumentParser(description='merge file statistics with file frequence generate by plink')
    parser.add_argument('--maf',type=float,default=0.01, help="basename of file for female")
    parser.add_argument('--base',type=str, help="basename of file for female")
    parser.add_argument('--miss',type=float,default=0.02, help="basename of file for female")
    parser.add_argument('--out', type=str,help="",default="out")
    args = parser.parse_args()
    return args
 
args = parseArguments()
 
def get_list_af(basename, maf) :
  data = pd.read_csv(basename+'.frq',delim_whitespace=True)
  data['balise_maf']=(data.loc[:,'MAF'].isna()==False) & (data.loc[:,'MAF'] >= maf) & (data.loc[:,'MAF'] <= (1-maf))
  return data

def get_missing(basename,rate) :
  data = pd.read_csv(basename+'.lmiss',delim_whitespace=True)
  data['balise_miss']=data['F_MISS']  < rate
  return data

def get_info_sex(basename, maf, ratemiss, sex) :
  data_af =get_list_af(basename, maf)
  data_miss =  get_missing(basename,ratemiss)
  dataall=data_af.merge(data_miss, suffixes=["_af", "_miss"], on=["CHR","SNP"])
  return dataall

 
data_all=get_info_sex(args.base, args.maf, args.miss, 'M')

data_all['N_A1']=data_all['MAF']*data_all['N_GENO']
data_all['N_A2']=data_all['N_GENO'] - data_all['N_A1']


data_all['all_filter'] =  True #data_all.balise_maf_male & data_all.balise_maf_female &
for head in ['balise_maf',  'balise_miss'] :
  data_all['all_filter']=(data_all[head] & data_all['all_filter'])
  
data_all.to_csv(args.out+'_resume.csv')

data_all.loc[data_all.all_filter,["SNP"]].to_csv(args.out+'.in',header=False, index=False)

pdfdic={}
pdfdic["ni"]=len(data_all.CHR)
pdfdic["nmal"]=data_all.NCHROBS.max()

pdfdic["filt_mal_maf"]=(data_all.balise_maf==False).sum()
pdfdic["filt_mal_miss"]=(data_all.balise_miss==False).sum()
pdfdic["allfilter"]=(data_all.all_filter==False).sum()

pdfdic["mafmale"] = args.maf
pdfdic["missmale"] = args.miss
pdfdic["resumecsv"]='*-protect*-url{'+args.out+'_resume.csv}'


template="""
Filters to analyse X chromosome are :
*-begin{itemize}
*-item individuals man are split
*-item Heterozygosity in males had been considered as missing.
*-end{itemize}



Data description :
*-begin{itemize}
*-item SNPs number on Y chromome (24) before filters : %(ni)s
*-item Male number : %(nmal)s
*-end{itemize}

*-begin{table}[h!]
*-centering
*-begin{tabular}{c|c}
 *-hline
 Filters  &  Snps number discarded *-*-
*-hline*-hline
 MAF *-textless %(mafmale)s in males & %(filt_mal_maf)s  *-*-
 Missingness *-textgreater  %(missmale)s in males & %(filt_mal_miss)s  *-*-
 Final filter & %(allfilter)s *-*-
*-hline
*-end{tabular}
*-caption{Filters resume on Y chromosome SNPs numbers deleted by filter }
*-end{table}

descriptive files of statistic for each position can be find in %(resumecsv)s
"""

letter=template%pdfdic

out=open("%s.tex"%args.out,"w")
out.write(letter)
out.close()



#datafemale_af=get_list_af(args.base_female, args.maf_female)
#datamale_af=get_list_af(args.base_male, args.maf_male)

# merge   CHR               SNP 




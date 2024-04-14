#!/usr/bin/env python3
#Load HWE P-value file and generate frequency_distribution
 
import pandas as   pd
import numpy  as np
import sys
from matplotlib import use
use('Agg')
import argparse
import matplotlib.pyplot as plt
import sys
 
def parseArguments():
    parser = argparse.ArgumentParser(description='merge file statistics with file frequence generate by plink')
    parser.add_argument('--base_female',type=str,required=True, help="basename of file for female")
    parser.add_argument('--maf_female',type=float,default=0.01, help="basename of file for female")
    parser.add_argument('--maf_male',type=float,default=0.01, help="basename of file for female")
    parser.add_argument('--miss_female',type=float,default=0.02, help="basename of file for female")
    parser.add_argument('--miss_male',type=float,default=0.02, help="basename of file for female")
    parser.add_argument('--diff_miss',type=float,default=0.02, help="basename of file for female")
    parser.add_argument('--base_male',type=str,required=True, help="basename of file for female")
    parser.add_argument('--threshold_hwe',type=float, default=1e-4, help="basename of file for female")
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

def get_hwe(basename, treshold) :
  data = pd.read_csv(basename+'.hwe',delim_whitespace=True)
  # CHR              SNP     TEST   A1   A2                 GENO   O(HET)   E(HET)            P 
  data['balise_hwe']=data['P']  >=treshold
  return data

def get_info_sex(basename, maf, ratemiss, sex) :
  data_af =get_list_af(basename, maf)
  data_miss =  get_missing(basename,ratemiss)
  dataall=data_af.merge(data_miss, suffixes=["_af", "_miss"], on=["CHR","SNP"])
  if sex == 'F' :
    data_hwe = get_hwe(basename, args.threshold_hwe)
    dataall=dataall.merge(data_hwe,  suffixes=["", "_hwe"],on=["CHR","SNP"])
  return dataall
  
 
datafemale=get_info_sex(args.base_female, args.maf_female, args.miss_female, 'F')
datamale=get_info_sex(args.base_male, args.maf_male, args.miss_male, 'M')

data_all=datafemale.merge(datamale,suffixes=["_female", "_male"], on=["CHR","SNP"])

## balise_miss_femal
data_all['balise_diff_miss']=(data_all.F_MISS_female - data_all.F_MISS_male).abs() < args.diff_miss
data_all['all_filter'] =  True #data_all.balise_maf_male & data_all.balise_maf_female &
for head in ['balise_maf_male', 'balise_maf_female', 'balise_miss_male', 'balise_miss_female', 'balise_hwe', 'balise_diff_miss'] :
  data_all['all_filter']=(data_all[head] & data_all['all_filter'])


data_all.to_csv(args.out+'_resume.csv')

data_all[["SNP"]].to_csv(args.out+'.in',header=False, index=False)

pdfdic={}
pdfdic["ni"]=len(data_all.CHR)
pdfdic["nmal"]=data_all.NCHROBS_male.max()
pdfdic["nfem"]=data_all.NCHROBS_female.max()/2

pdfdic["filt_mal_maf"]=(data_all.balise_maf_male==False).sum()
pdfdic["filt_femal_maf"]=(data_all.balise_maf_female==False).sum()
pdfdic["filt_mal_miss"]=(data_all.balise_miss_male==False).sum()
pdfdic["filt_femal_miss"]=(data_all.balise_miss_female==False).sum()
pdfdic["filt_femal_hwe"]=(data_all.balise_hwe==False).sum()
pdfdic["filt_diff_miss"]=(data_all.balise_diff_miss==False).sum()
pdfdic["allfilter"]=(data_all.all_filter==False).sum()

pdfdic["mafmale"] = args.maf_male
pdfdic["maffemale"] = args.maf_female
pdfdic["missmale"] = args.miss_male
pdfdic["missfemale"] = args.miss_female
pdfdic["threshold_hwe"] = args.threshold_hwe
pdfdic["diffmiss"] = args.diff_miss
pdfdic["resumecsv"]='*-protect*-url{'+args.out+'_resume.csv}'



template="""
Filters to analyse X chromosome are :
*-begin{itemize}
*-item based on papers from Inke R. KOnig, Christina Loley, Jeanette Erdmann, Andreas Ziegler,  How to Include Chromosome X in Your Genome-Wide Association Study
*-item option `--split-x` used to split chromosome X and pseudo X
*-item individuals are split by sex 
*-item Heterozygosity in males had been considered as missing.
*-end{itemize}



Data description :
*-begin{itemize}
*-item SNPs number on X chromome (26) before filters : %(ni)s
*-item Female number : %(nfem)s
*-item Male number : %(nmal)s
*-end{itemize}

*-begin{table}[h!]
*-centering
*-begin{tabular}{ c | c }
 *-hline
 Filters  &  Snps number discarded *-*-
*-hline*-hline
 MAF *-textless %(mafmale)s in males & %(filt_mal_maf)s  *-*-
 MAF *-textless %(maffemale)s in females & %(filt_femal_maf)s *-*- 
 Missingness *-textgreater  %(missmale)s in males & %(filt_mal_miss)s  *-*-
 Missingness *-textgreater %(missfemale)s in females & %(filt_femal_miss)s *-*- 
 P HWE *-textless %(threshold_hwe)s in females & %(filt_femal_hwe)s  *-*-
 abs(Mif Male - Mif female) *-textgreater %(diffmiss)s in females & %(filt_diff_miss)s *-*-
 Final filter & %(allfilter)s *-*-
*-hline
*-end{tabular}
*-caption{Filters resume on X chromosome SNPs numbers deleted by filter }
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




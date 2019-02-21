#!/usr/bin/env python3

# (c) University of the Witwatersand, Johannesburg on behalf of the H3ABioinformatics Network Consortium
# 2016-2018
# Licensed under the Creative Commons Attribution 4.0 International Licence. 
# See the "LICENSE" file for details

import pandas as pd
import sys
import numpy as np
import glob
import matplotlib as mpl
mpl.use('Agg') 
import matplotlib.pyplot as plt
import os
import re


g = glob.glob


label   = sys.argv[1]
base    = sys.argv[2]
test    = sys.argv[3]
pheno   = sys.argv[4]
gotcovar = sys.argv[5]
gtype    = sys.argv[6]

if len(pheno)==0: pheno = "infam"

EOL = chr(10)

def drawManhatten(pheno, result,outf):
    fig, (ax1, ax2) =  plt.subplots(2, 1, sharey=True)
    plt.xlabel("Chromosome",fontsize=14)
    for ax in [ax1,ax2]:
        ax.set_ylabel("$-\log P$",fontsize=14)
        ax.xaxis.set_tick_params(labelsize=11)
        ax.xaxis.set_tick_params(labelsize=11)
        ax.yaxis.set_tick_params(labelsize=11)
        ax.yaxis.set_tick_params(labelsize=11)
    ax=ax1
    delta=0
    colours= ['crimson','blue','green']
    xtick_pos=[]
    xtick_label = []
    chroms = result.groupby("CHR")
    for chrom_num, chrom_res in chroms:
        this_chrom = result['CHR']==chrom_num
        result.loc[this_chrom,'BP']+=delta
        old_delta=delta
        delta = delta + int(chrom_res.tail(1)['BP'])
        xtick_pos.append((delta+old_delta)/2)
        xtick_label.append(str(chrom_num))
        under_thresh = result['P']<0.005
        ax.scatter(result.loc[this_chrom & under_thresh, 'BP'],\
                   -np.log10(result.loc[this_chrom  & under_thresh,'P']),c=colours[chrom_num%3])
        if chrom_num == 9:
           ax.set_xticklabels(xtick_label)
           ax.set_xticks(xtick_pos)
           xtick_pos=[]
           xtick_label=[]
           ax=ax2
           delta=0
    ax.set_xticklabels(xtick_label)
    ax.set_xticks(xtick_pos)
    plt.savefig(outf)

def drawQQ(pheno,result,outf):
    plt.figure()
    sort_p = -np.log10(result['P'].sort_values())
    n=len(sort_p)
    expected = -np.log10(np.linspace(1/n,1,n))
    plt.plot(expected,sort_p)
    plt.xlabel("Expected $-\log P$-value",fontsize=14)
    plt.ylabel("Observed $-\log P$-value",fontsize=14)
    plt.plot(expected,expected)
    plt.savefig(outf)


qq_template = """


*-section{PLINK Results for phenotype *-protect*-url{%(pheno)s}, test %(test)s}

The raw results can be found in the *-textbf{%(test)s} directory.
The QQ plot can be found in Figure *-ref{fig:qq:%(pheno)s:%(test)s}, and the Manhatten
plot in Figure *-ref{fig:man:%(pheno)s:%(test)s}.

"""


qq_figs = """
*-begin{figure}[hb]
*-begin{center}
*-includegraphics[width=17cm]{%(qqfile)s}
*-end{center}
*-caption{QQ plot for PLINK testing --  *-protect*-url{%(pheno)s}, test %(test)s -- *-protect*-url{%(testing)s} }
*-label{fig:qq:%(pheno)s:%(test)s}
*-end{figure}

*-begin{figure}[ht]
*-begin{center}
*-includegraphics[width=17cm]{%(manfile)s}
*-end{center}
*-caption{PLINK testing: Manhatten plot for --  *-protect*-url{%(pheno)s}, test %(test)s -- *-protect*-url{%(testing)s}}
*-label{fig:man:%(pheno)s:%(test)s}
*-end{figure}
"""

interesting=r"""


*-item The top ranking SNPs for %s according to %s are found in Table~\ref{tab:plink:top:%s}

*-begin{table}
*-begin{center}
*-begin{tabular}{l l}*-hline
SNP & P *-*-*-hline
"""

endtab = r"""
*-end{tabular}
*-end{center}
*-caption{Top SNPs found by plink using %s for %s test *-protect*-url{%s}}
*-label{tab:plink:top:%s}
*-end{table}
"""

def lambdaGC(log):
    if not os.path.exists(log): return ""
    res = ""
    f = open(log)
    for line in f:
        m = re.search("Genomic inflation est. lambda(.*)", line)
        if m:
            res = "*-paragraph*{Genomic inflation: }: The estimate of genomic inflation ##*-lambda## "+m.group(1)+EOL+EOL
            break
    return res


def showInteresting(frame, p_col, test_name):
    cpheno = pheno.replace("_","-")
    slist = interesting % (str(p_col),test_name,str(p_col)+cpheno)
    for row in frame.iterrows():
        slist = slist+"%s & %6.4f "%(row[1]['SNP'].replace("_","-"),row[1][p_col])+r"*-*-"+chr(10)
    slist = slist + endtab%(test_name,p_col,pheno,str(p_col)+cpheno)
    return slist

def processProbs(data,p_col,p_name):
    mf = pd.read_csv(data,delim_whitespace=True)
    mf.sort_values(by=p_col,inplace=True)
    top10= mf.head(n=10)
    selected = mf[mf[p_col]<=0.05]
    if len(selected)>10:
        selected = top10
    slist = showInteresting(selected,p_col,p_name)
    return slist


def clean(name):
   return name.replace("/","-").replace("np.","")

def showResults():
    # first get the result for the straight forward test
    if test=="assoc":
        data = clean(pheno) + ".*assoc"
    else:
        data = clean(pheno) + ".*assoc."+test
    data = glob.glob(data)
    if len(data)!=1:
        print((EOL*3)+"---- Expected only one assoc file found <"+",".join(result)+">"+(EOL*3))
        sys.exit(17)
    data = data[0]
    asf   = pd.read_csv(data,delim_whitespace=True)
    if gotcovar == "1" and test != "assoc":
        asf = asf[asf['TEST']=="ADD"]
    cpheno  = pheno.replace("_","-")
    hashd = {'manfile':clean("%s-pl-man-%s.%s"%(base,cpheno,gtype)), 'testing':data,'test':test, 'pheno':cpheno,
             'qqfile':clean("%s-pl-qq-%s.%s"%(base,cpheno,gtype))}
    drawManhatten(pheno.replace("_","-"),asf,hashd['manfile'])
    drawQQ(pheno.replace("_","-"),asf,hashd['qqfile'])
    outpics = qq_template%hashd + lambdaGC(clean(pheno)+".log")+EOL+EOL+"*-begin{itemize}"+EOL
    for (correct,col,label) in [(".mperm",'EMP2',"permutation testing"),(".adjusted",'BONF',"Bonferroni correction")]:
        nd = data + correct
        if os.path.exists(nd):
            outpics = outpics+processProbs(nd,col,label)
    outpics =outpics+"*-end{itemize}"+EOL+EOL+qq_figs%hashd+EOL+"*-clearpage"+EOL
    return outpics
        
out_pics = showResults()



textf = open("%s%s%s.tex"%(label,pheno.replace("_","-"),test),"w")
textf.write(out_pics.replace("*-",chr(92)).replace("##",chr(36)).replace("@.@",chr(10)))
textf.close()

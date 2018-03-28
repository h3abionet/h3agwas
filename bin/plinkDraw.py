#!/usr/bin/env python3

import pandas as pd
import sys
import numpy as np
import glob
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import os
g = glob.glob


base    = sys.argv[1]
test    = sys.argv[2]
phenos  = sys.argv[3].split(",")
gotcovar = sys.argv[4]
gtype    = sys.argv[5]

if len(phenos[0])==0:
    psep = ""
else:
    psep = "."

EOL = chr(10)

def drawManhatten(pheno, result,outf):
    fig, (ax1, ax2) =  plt.subplots(2, 1, sharey=True)
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
    plt.plot(expected,expected)
    plt.savefig(outf)


qq_template = """

\section{PLINK Results for phenotype *-protect*-url{%(pheno)s}, test %(test)s}

The QQ plot can be found in Figure *-ref{fig:qq}, and the Manhatten
plot in Figure *-ref{fig:man}.

*-begin{figure}[ht]
*-begin{center}
*-includegraphics[width=10cm]{%(qqfile)s}
*-end{center}
*-caption{QQ plot for PLINK testing --  *-protect*-url{%(pheno)s}, test %(test)s -- *-protect*-url{%(testing)s} }
*-label{fig:qq}
*-end{figure}

*-begin{figure}[ht]
*-begin{center}
*-includegraphics[width=10cm]{%(manfile)s}
*-end{center}
*-caption{PLINK testing: Manhatten plot for --  *-protect*-url{%(pheno)s}, test %(test)s -- *-protect*-url{%(testing)s}}
*-label{fig:man}
*-end{figure}
"""

interesting=r"""

The top ranking SNPs according to %s are found in Table~\ref{tab:plink:top:%s}

*-begin{table}
*-begin{center}
*-begin{tabular}{l l}*-hline
SNP & P *-*-*-hline
"""

endtab = r"""
*-end{tabular}
*-end{center}
*-caption{Top SNPs found by plink using %s for %s}
*-label{tab:plink:top:%s}
*-end{table}
"""

def showInteresting(frame, p_col, test_name):
    slist = interesting % ("permutation testing",str(p_col)+test_name)
    for row in frame.iterrows():
        slist = slist+"%s & %6.4f "%(row[1]['SNP'].replace("_","-"),row[1]['EMP2'])+r"*-*-"+chr(10)
    slist = slist + endtab%(test_name,p_col,str(p_col)+test_name)
    return slist

def processProbs(data,p_col,p_name):
    mf = pd.read_csv(data,delim_whitespace=True)
    mf.sort_values(by=p_col,inplace=True)
    top10= mf.head()
    selected = mf[mf[p_col]<=0.05]
    if len(selected)<10:
        selected = top10
    slist = showInteresting(selected,p_col,p_name)
    return slist

def showResults():
    outpics = ""
    for pheno in phenos:
        # first get the result for the straight forward test
        data = base+psep + pheno.split("/")[0]+".assoc."+test
        if len(pheno)==0: pheno = "in fam"
        asf   = pd.read_csv(data,delim_whitespace=True)
        if gotcovar == "1":
            asf = asf[asf['TEST']=="ADD"]
        hashd = {'manfile':("%s-man-%s.%s"%(base,pheno,gtype)).replace("/","-").replace("np.",""), 'testing':data,'test':test, 'pheno':pheno,\
                 'qqfile':("%s-qq-%s.%s"%(base,pheno,gtype)).replace("/","-").replace("np.","")}
        drawManhatten(pheno,asf,hashd['manfile'])
        drawQQ(pheno,asf,hashd['qqfile'])
        outpics = outpics  + qq_template%hashd
        for (correct,col,label) in [(".mperm",'EMP2',"permutation testing"),(".adjusted",'BONF',"Bonferroni correction")]:
            nd = data + correct
            if os.path.exists(nd):
                outpics = outpics+processProbs(nd,col,label)
        return outpics
        
out_pics = showResults()



textf = open("C050%s.tex"%test,"w")
textf.write(out_pics.replace("*-",chr(92)).replace("##",chr(36)).replace("@.@",chr(10)))
textf.close()

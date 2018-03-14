#!/usr/bin/env python3

import pandas as pd
import sys
import numpy as np
import glob
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt


assoc = glob.glob("*assoc")
mperm = glob.glob("*mperm")
adjt  = glob.glob("*adjusted")

out_man = sys.argv[1]
out_qq  = sys.argv[2]


def drawManhatten(result):
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
    plt.savefig(out_man)

def drawQQ(result):
    plt.figure()
    sort_p = -np.log10(result['P'].sort_values())
    n=len(sort_p)
    expected = -np.log10(np.linspace(1/n,1,n))
    plt.plot(expected,sort_p)
    plt.plot(expected,expected)
    plt.savefig(out_qq)


qq_template = """

The QQ plot can be 
found in Figure *-ref{fig:qq}, and  the Manhatten plot in Figure *-ref{fig:man}.

*-begin{figure}[ht]
*-includegraphics[width=10cm]{%s}
*-caption{QQ plot}
*-label{fig:qq}
*-end{figure}

*-begin{figure}[ht]
*-includegraphics[width=10cm]{%s}
*-caption{Manhatten plot}
*-label{fig:man}
*-end{figure}
"""

interesting=r"""

The top ranking SNPs according to %s are
*-begin{tabular}{l l}*-hline
SNP & P *-*-*-hline
"""


def showInteresting(frame, p_col, test_name):
    slist = interesting % "permutation testing"
    for row in frame.iterrows():
        slist = slist+"%s & %6.4f "%(row[1]['SNP'],row[1]['EMP2'])+r"*-*-"+chr(10)
    slist = slist + r"*-end{tabular}"
    return slist

def processProbs(frame,p_col,p_name):
    mf = pd.read_csv(mperm[0],delim_whitespace=True)
    mf.sort_values(by=p_col,inplace=True)
    top10= mf.head()
    selected = mf[mf[p_col]<=0.05]
    if len(selected)<10:
        selected = top10
    slist = showInteresting(selected,p_col,p_name)
    return slist

if len(assoc)> 0:
    asf   = pd.read_csv(assoc[0],delim_whitespace=True)
    drawManhatten(asf)
    drawQQ(asf)
    out_pics = qq_template%(out_qq,out_man)
else:
    out_pics = ""

slist =""
if len(mperm)>0:
    slist = processProbs(mperm[0],'EMP2',"permutation testing")
if len(adjt)>0:
    slist = slist + processProbs(mperm[0],'BONF',"Bonferroni correction")

textf = open(sys.argv[3],"w")
textf.write(out_pics.replace("*-",chr(92)).replace("##",chr(36)).replace("@.@",chr(10)))
textf.write(slist.replace("*-",chr(92)).replace("##",chr(36)).replace("@.@",chr(10)))
textf.close()

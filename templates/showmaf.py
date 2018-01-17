#!/usr/bin/env python3

import pandas as pd
import numpy  as np
import sys
from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt

EOL=chr(10)

if len(sys.argv)<=1:
    sys.argv="showmaf.py $freq  $base".split()



template = """
*-paragraph*{Minor allele frequency.} Table *-ref{tab:initmafspec} on page *-pageref{tab:initmafspec} shows the minor allele frequency spectrum for the raw data. The number of monomorphic SNPs is shown in the first row. Note that some of the MAFs with very low MAF are actually monomorphic, with the polymorphisms due to genotyping error. Figure *-ref{fig:initmafspec} on page *-pageref{fig:initmafspec} shows the cumulative distribution of MAF.  This can be used to determine an appropriate MAF cut-off.


Note that the *-emph{minor} allele is determined with respect to the frequency spectrum inthis data -- `minor' is not synonym for alternate or non-reference allele, or the allele that has minor frequency in some other data set. Under this definition the MAF is always ##*-leq 0.5##.


*-begin{table}[htb]*-centering
%s
*-caption{Minor Allele Frequency spectrum of the raw data. The number of apparently monommorphic SNPs is shown in the row labelled 0; the other rows show the number of SNPs in the bins shown.}
*-label{tab:initmafspec}
*-end{table}

*-begin{figure}[htb]*-centering
*-includegraphics[width=10cm]{%s}
*-caption{Cumulative frequency of SNPs. For a frequency shown on the ##x##-axis, the corresponding ##y##-value shows the proportion of SNPs with frequency *-emph{at least this frequency}; that is, it shows the proportion of SNPs which will *-emph{remain} if the MAF filter of this ##x##-value is chosne.}
*-label{fig:initmafspec}
*-end{figure}
"""

def rename(x):
    if x.left < 0:
       return "0"
    else:
       return str(x)  


def getTable(frm):
   frm["MAF bin"]=pd.cut(frm['MAF'],[-0.001,0,0.005,0.01,0.02,0.03,0.04,0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50],right=True)
   g1 = frm.groupby('MAF bin').size()
   g1 = g1.rename(rename)
   g1.name = "Num SNPs"
   return g1.to_latex(bold_rows=True)


def getPic(frm,fname):
    r = [-0.00001]+list(map(lambda x: x/100,range(0,51,2)))
    xs = list(map(lambda x:x/100,range(0,51,2)))
    frm['bin']=pd.cut(frm['MAF'],r,right=True)
    g2  = frm.groupby('bin').size()
    cum = (g2.sum()-g2.cumsum()+g2.iloc[0])/g2.sum()
    plt.plot(xs,cum)
    plt.xlim(0,0.5)
    plt.ylabel('Proportion SNPs ##*-geq## this freq'.replace("*-",chr(92)).replace("##",chr(36)))
    plt.xlabel('Frequency')
    plt.savefig(fname)



frm = pd.read_csv(sys.argv[1],delim_whitespace=True)
base = sys.argv[2]
pdfout = "%s.pdf"%base
texout = "%s.tex"%base

getPic(frm,pdfout)
curr = getTable(frm)


g=open(texout,"w")
g.write((template%(curr,pdfout)).replace("*-",chr(92)).replace("##",chr(36)))
g.close()

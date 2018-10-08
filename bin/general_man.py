#!/usr/bin/env python3

''' 
Produces Manhatten, QQ plot and supporting tex file  for output from some tools
'''


import argparse
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys


def parseArguments():
    parser = argparse.ArgumentParser(description='Produces Manhatten, QQ plot and supporting tex file  for output from some tools')
    parser.add_argument('--inp',type=str,required=True, help="association files")
    parser.add_argument('--phenoname',type=str,required=True,help="pheno names")
    parser.add_argument('--out', type=str,help="out of tex file",default="")
    parser.add_argument('--chro_header',type=str,required=True,help="chro header in inp files")
    parser.add_argument('--pos_header',type=str,required=True,help="pos header in inp files")
    parser.add_argument('--rs_header',type=str,required=True,help="rs header in inp files")
    parser.add_argument('--pval_header',type=str,required=True,help="pvalue header in inp files")
    parser.add_argument('--beta_header',type=str,required=True,help="beta header in inp files")
    parser.add_argument('--info_prog',type=str,required=True,help="info program to be print in tex files")
    args = parser.parse_args()
    return args

args = parseArguments()

inp     = args.inp
phenos= args.phenoname
base    = args.out.replace("/np.","-").replace("/","").replace(".","-").replace("_","-")
ChroEnt=args.chro_header
PosEnt=args.pos_header
RsEnt=args.rs_header
PvalueEnt=args.pval_header
BEnt=args.beta_header
NameProg=args.info_prog

out_man = "%s-man.pdf"%base
out_qq  = "%s-qq.pdf"%base
out_tex = "%s.tex"%base


EOL = chr(10)
latex_tex = r"""
\clearpage
\section{Result of %(NameProg)s analysis : phenotype %(pheno)s}
\label{sec:%(NameProg)s}

All the results from the %(NameProg)s analysis can be found in the \textbf{%(NameProg)s} directory.
The result of the %(NameProg)s analysis is shown for phentoytype \protect\url{%(pheno)s}.  The file with association statistics is found in \protect\url{%(fname)s}. The top 10 results are shown in Table~\ref{tab:top:%(pheno)s}:

\begin{table}[tb]
\begin{center}
\begin{tabular}{|l l l l r|}\hline
Chr & SNP & Pos & Beta & P \\\hline
%(best)s\hline
\end{tabular}
\end{center}
\caption{The top 10 SNPs found by %(NameProg)s analysis of phenotype \protect\url{%(pheno)s}}
\label{tab:top:%(pheno)s}
\end{table}

\noindent
The Manhatten plot can be found in Figure \ref{fig:%(NameProg)s-man-%(pheno)s}. The corresponding QQ-plot can
be found in Figure~\ref{fig:%(NameProg)s-qq-%(pheno)s}. 

\begin{figure}[ht]
\begin{center}
\includegraphics[width=15cm]{%(man)s}
\caption{%(NameProg)s testing: Manhatten plot for phenotype %(pheno)s}
\label{fig:%(NameProg)s-man-%(pheno)s}
\end{center}
\end{figure}

\begin{figure}[ht]
\begin{center}
\includegraphics[width=15cm]{%(qq)s}
\caption{%(NameProg)s testing: QQ-plot for phenotype %(pheno)s}
\label{fig:%(NameProg)s-qq-%(pheno)s}
\end{center}
\end{figure}

\clearpage
"""

best_row = r"%s & %s & %s & %7.4f & %6.3E \\"

def get10Best(result,pheno):
    top = result.sort_values(PvalueEnt).head(10)
    best = ""
    for r in top.iterrows():
        x= r[1]
        row = best_row%(x[ChroEnt],x[RsEnt].replace("_","-"),x[PosEnt],x[BEnt],x[PvalueEnt])+EOL
        best = best + row
    return best
        

result = pd.read_csv(inp,delim_whitespace=True)
chroms = result.groupby(ChroEnt)

fig, (ax1, ax2) =  plt.subplots(2, 1, sharey=True)

ax=ax1
delta=0
colours= ['crimson','blue','green']
xtick_pos=[]
xtick_label = []
result['ps_new']=result[PosEnt]
for chrom_num, chrom_res in chroms:
    this_chrom = result[ChroEnt]==chrom_num
    result.loc[this_chrom,'ps_new']+=delta
    old_delta=delta
    delta = delta + int(chrom_res.tail(1)['ps_new'])
    xtick_pos.append((delta+old_delta)/2)
    xtick_label.append(str(chrom_num))
    under_thresh = result[PvalueEnt]<0.005
    ax.scatter(result.loc[this_chrom & under_thresh, 'ps_new'],\
               -np.log10(result.loc[this_chrom  & under_thresh,PvalueEnt]),c=colours[chrom_num%3])
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

plt.figure()
sort_p = -np.log10(result[PvalueEnt].sort_values())
n=len(sort_p)
if n==0 :
   expected=[]
else :
   expected = -np.log10(np.linspace(1/n,1,n))
plt.plot(expected,sort_p)
plt.plot(expected,expected)
plt.savefig(out_qq)



tex_file=""
#pheno=phenos[1]
pheno=phenos
best = get10Best(result,pheno)
hashd = { 'pheno':pheno.replace("_","-"), 'man':out_man, 'qq':out_qq, 'best':best, 'fname':inp, 'NameProg':NameProg}
tex_file = tex_file + latex_tex%hashd

g = open(out_tex,"w")
g.write(tex_file)
g.close()



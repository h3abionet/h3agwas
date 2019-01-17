#!/usr/bin/env python3


import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

inp     = sys.argv[1]
phenos  = sys.argv[2].split("-")
base    = sys.argv[3].replace("/np.","-").replace("/","").replace(".","-").replace("_","-")
out_man = "%s-gemma-man.pdf"%base
out_qq  = "%s-gemma-qq.pdf"%base
out_tex = "%s.tex"%base


EOL = chr(10)
C049 = r"""
\clearpage
\section{Result of Gemma analysis : phenotype %(pheno)s}
\label{sec:gemma}

All the results from the GEMMA analysis can be found in the \textbf{gemma} directory.
The result of the GEMMA analysis is shown for phentoytype \protect\url{%(pheno)s}.  The file with association statistics is found in \protect\url{%(fname)s}. The top 10 results are shown in Table~\ref{tab:top:%(pheno)s}:

\begin{table}[tb]
\begin{center}
\begin{tabular}{|l l l l r|}\hline
Chr & SNP & Pos & Beta & P \\\hline
%(best)s\hline
\end{tabular}
\end{center}
\caption{The top 10 SNPs found by GEMMA analysis of phenotype \protect\url{%(pheno)s}}
\label{tab:top:%(pheno)s}
\end{table}

\noindent
The Manhatten plot can be found in Figure \ref{fig:gemma-man-%(pheno)s}. The corresponding QQ-plot can
be found in Figure~\ref{fig:gemma-qq-%(pheno)s}. 

\begin{figure}[ht]
\begin{center}
\includegraphics[width=15cm]{%(man)s}
\caption{Gemma testing: Manhatten plot for phenotype %(pheno)s}
\label{fig:gemma-man-%(pheno)s}
\end{center}
\end{figure}

\begin{figure}[ht]
\begin{center}
\includegraphics[width=15cm]{%(qq)s}
\caption{Gemma testing: QQ-plot for phenotype %(pheno)s}
\label{fig:gemma-qq-%(pheno)s}
\end{center}
\end{figure}

\clearpage
"""

best_row = r"%s & %s & %s & %7.4f & %6.3E \\"

def get10Best(result,pheno):
    top = result.sort_values('p_wald').head(10)
    best = ""
    for r in top.iterrows():
        x= r[1]
        row = best_row%(x['chr'],x['rs'].replace("_","-"),x['ps'],x['beta'],x['p_wald'])+EOL
        best = best + row
    return best
        

result = pd.read_csv(inp,delim_whitespace=True)
chroms = result.groupby("chr")

fig, (ax1, ax2) =  plt.subplots(2, 1, sharey=True)

ax=ax1
delta=0
colours= ['crimson','blue','green']
xtick_pos=[]
xtick_label = []
result['ps_new']=result['ps']
for chrom_num, chrom_res in chroms:
    this_chrom = result['chr']==chrom_num
    result.loc[this_chrom,'ps_new']+=delta
    old_delta=delta
    delta = delta + int(chrom_res.tail(1)['ps_new'])
    xtick_pos.append((delta+old_delta)/2)
    xtick_label.append(str(chrom_num))
    under_thresh = result['p_wald']<0.005
    ax.scatter(result.loc[this_chrom & under_thresh, 'ps_new'],\
               -np.log10(result.loc[this_chrom  & under_thresh,'p_wald']),c=colours[chrom_num%3])
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
sort_p = -np.log10(result['p_wald'].sort_values())
n=len(sort_p)
expected = -np.log10(np.linspace(1/n,1,n))
plt.plot(expected,sort_p)
plt.plot(expected,expected)
plt.savefig(out_qq)



tex_file=""
pheno=phenos[1]
best = get10Best(result,pheno)
hashd = { 'pheno':pheno.replace("_","-"), 'man':out_man, 'qq':out_qq, 'best':best, 'fname':inp}
tex_file = tex_file + C049%hashd

g = open(out_tex,"w")
g.write(tex_file)
g.close()



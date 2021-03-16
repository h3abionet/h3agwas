#!/usr/bin/env python3

''' 
Produces  QQ plot and supporting tex file  for output from some tools
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
    parser.add_argument('--out', type=str,help="out of tex file",default="")
    parser.add_argument('--rs_header',type=str,required=True,help="rs header in inp files")
    parser.add_argument('--pval_header',type=str,required=True,help="pvalue header in inp files")
    parser.add_argument('--beta_header',type=str,required=True,help="beta header in inp files")
    parser.add_argument('--info_prog',type=str,required=True,help="info program to be print in tex files")
    args = parser.parse_args()
    return args

args = parseArguments()

inp     = args.inp
out     = args.out
RsEnt=args.rs_header
PvalueEnt=args.pval_header
BEnt=args.beta_header
NameProg=args.info_prog

out_qq  = "%s-qq.pdf"%out
out_tex = "%s.tex"%out


EOL = chr(10)
latex_tex = r"""
\clearpage
\section{Result of meta analysis : %(NameProg)s}
\label{sec:%(NameProg)s}

All the results from the %(NameProg)s analysis can be found in the \textbf{%(NameProg)s} directory.
The result of the %(NameProg)s analysis is shown for program \protect\url{%(NameProg)s}.  The file with association statistics is found in \protect\url{%(fname)s}. The top 10 results are shown in Table~\ref{tab:top:%(NameProg)s}:

\begin{table}[tb]
\begin{center}
\begin{tabular}{|l l l l r|}\hline
SNP & Beta & P \\\hline
%(best)s\hline
\end{tabular}
\end{center}
\caption{The top 10 SNPs found by %(NameProg)s analysis}
\label{tab:top:%(NameProg)s}
\end{table}

\noindent

\begin{figure}[ht]
\begin{center}
\includegraphics[width=15cm]{%(qq)s}
\caption{%(NameProg)s testing: QQ-plot for %(NameProg)s}
\label{fig:%(NameProg)s-qq}
\end{center}
\end{figure}

\clearpage
"""

best_row = r"%s & %7.4f & %6.3E \\"

def get10Best(result):
    top = result.sort_values(PvalueEnt).head(10)
    best = ""
    for r in top.iterrows():
        x= r[1]
        row = best_row%(x[RsEnt],x[BEnt],x[PvalueEnt])+EOL
        best = best + row
    return best

result = pd.read_csv(inp,delim_whitespace=True, keep_default_na=True, na_values=["nan","nane-nan", "NA"])

result[PvalueEnt] = result[PvalueEnt].astype(float)
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
best = get10Best(result)
hashd = {'qq':out_qq, 'best':best, 'fname':inp, 'NameProg':NameProg}
tex_file = tex_file + latex_tex%hashd

g = open(out_tex,"w")
g.write(tex_file)
g.close()



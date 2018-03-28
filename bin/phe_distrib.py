#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy  as np
import sys
import matplotlib
from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt

EOL=chr(10)


def parseArguments():
    if len(sys.argv)<=1:
        sys.argv="mafplot.py $input $output".split()
    parser=argparse.ArgumentParser()
    parser.add_argument("--phenos", type=str, metavar='phenotypes', required=True)
    parser.add_argument("--skip-zero",  dest="skip_zero", action="store_true", default=False)
    parser.add_argument('input', type=str, metavar='input'),
    parser.add_argument('output', type=str, metavar='output'),
    args = parser.parse_args()
    return args





transforms = [[1,np.log1p],[np.sqrt,np.cbrt]]
transform_names = [["no transform","log transform"],["square root transform","cube root transform"]]


def numfrm(x):
    xstr=str(x)
    if "." not in xstr: return xstr
    if x<0.1:
        xstr="%6.4E"%x
    else:
        xstr = str(x)
        xstr = xstr[:xstr.index(".")+3]
    return xstr


def summary2LaTeX(summary,output,suf,pheno):
     phelab = pheno.replace("_",":")
     lat = EOL+EOL+\
     r"\begin{table}[hb]"+EOL+\
     r"\begin{center}"+EOL+r"\begin{tabular}{l  r  D{.}{.}{3}  D{.}{.}{3}   D{.}{.}{4}   D{.}{.}{4}} \\"+EOL + \
     r"Data & Count & \multicolumn{1}{c}{Min} & \multicolumn{1}{c}{Max} & \multicolumn{1}{c}{Ave} & \multicolumn{1}{c}{StdDev} \\\hline" +EOL
     for s in summary:
        lat = lat+" & ".join([s[0]]+list(map(numfrm,[s[1].count(),s[1].min(),s[1].max(),s[1].mean(),s[1].std()])))+ r"\\"+EOL
     lat = lat + r"\hline\end{tabular}"+EOL+r"\end{center}"+EOL+(r"""
       *-caption{Overview of phenotype *-protect*-url{%s} distribution}
       *-label{tab:overview:%s}
       *-end{table}
     """)%(pheno,phelab)
     lat = r"""

        A summary of the data for \url{%s} can be found in the Table~*-ref{tab:overview:%s}, transformed using 
        different transforms. A histogram is found in Figure \ref{fig:%s}.

        """ + lat + r"""

        
        \ourfig{fig:%s}{Histogram of *-protect*-url{%s} values under different transforms}{%s.%s}


        """
     return lat%(pheno,phelab,output,output,pheno,output,suf)

def showPheno(pname,frm):
    if args.skip_zero:
        data  = frm[frm[pname]>0][pname]
    else:
        data = frm[pname]
    fig,axs = plt.subplots(2,2)
    matplotlib.rcParams['xtick.labelsize']=13
    matplotlib.rcParams['ytick.labelsize']=13
    summary=[]
    for r in range(2):
        for c in range(2):
            axs[r][c].set_xlabel(transform_names[r][c],fontsize=12)
            axs[r][c].set_ylabel("Frequency",fontsize=12)
            fn = transforms[r][c]
            pdata = fn(data) if fn != 1 else data
            pdata = pdata[pdata.notnull()]
            summary.append((transform_names[r][c],pdata))
            axs[r][c].hist(pdata,bins=100)
    plt.tight_layout()
    output = ("%s-%s"%(args.output,pname)).replace("_","-")
    plt.savefig("%s.pdf"%output)
    return summary2LaTeX(summary,output,"pdf",pname)

args=parseArguments()
frm = pd.read_csv(args.input,delim_whitespace=True)
phenos = args.phenos.split(",")
output_latex= ""
for phen in phenos:
    dets = phen.split("/")
    pname = dets[0]
    output_latex = output_latex + showPheno(pname,frm)

g = open("%s.tex"%args.output,"w")
g.write(output_latex.replace("*-",chr(92)).replace("##",chr(36)))
g.close()




#!/usr/bin/env python3
''' 
Produces report for rs file
'''
from  subprocess import CalledProcessError
import subprocess
import os
import argparse
import numpy as np
import pandas as pd
import sys

EOL = chr(10)

def check_output(x,shell=False):
   ans=subprocess.check_output(x,shell=shell)
   ans=str(ans,'ascii') # encoding option only frmo python 3.6
   return ans


try:
   kpsewhich=check_output("which kpsewhich",shell=True)
   if kpsewhich:
      kpsewhich=check_output("kpsewhich datetime.sty",shell=True)
except CalledProcessError:
   kpsewhich=""


fancy="""
*-usepackage{fancyhdr}
*-usepackage[yyyymmdd,hhmmss]{datetime}
*-pagestyle{fancy}
*-rfoot{Completed on *-today*- at *-currenttime}
*-cfoot{}
*-lfoot{Page *-thepage}
"""

dateheader=""
if len(kpsewhich)>1:
   dfmt = kpsewhich.rstrip()
   if os.access(dfmt,os.R_OK):
      with open(dfmt) as f:
         for line in f:
            m=re.search("ProvidesPackage.datetime..(..../../..)",line)
            if m and m.group(1) >= "2010/09/21":
               dateheader=fancy



template='''
*-documentclass[11pt]{article}

*-usepackage[paper=a4paper,left=2cm,right=2cm,top=2cm,bottom=2cm]{geometry}
*-usepackage{graphicx}
*-usepackage{subfig}
*-usepackage{listings}
*-usepackage{longtable}
*-usepackage{array}
*-usepackage{booktabs}
*-usepackage{float}
*-usepackage{pdfpages}
*-usepackage{dcolumn}
*-floatstyle{ruled}
*-restylefloat{figure}
*-restylefloat{table}
*-newcommand{*-lefttblcol}{*-raggedright*-hspace{0pt}}
*-newcommand{*-righttblcol}{*-raggedleft*-hspace{0pt}}
*-newcommand{*-centretblcol}{*-centering*-hspace{0pt}}

*-newcolumntype{P}[1]{>{*-lefttblcol}p{#1}}
*-newcolumntype{Q}[1]{>{*-righttblcol}p{#1}}
*-newcolumntype{R}[1]{>{*-centretblcol}p{#1}}
*-lstset{
basicstyle=*-small*-ttfamily,
columns=flexible,
breaklines=true
}
*-usepackage{url}
*-title{Annotation Testing  %(rsname)s : %(pheno)s}
*-date{%(date)s}
'''+dateheader+'''
*-author{H3Agwas annotation Testing Pipeline}

*-begin{document}

*-maketitle

*-section{Introduction}

This report gives a brief overview of the run of the annotation of rs testing pipeline.
*-begin{itemize}
*-item You were testing for the following rs *-url{%(rsname)s}
*-item You were using the following pheno [*-url{%(pheno)s}]
*-item You were using the following covariates [*-url{%(cov)s}]
*-end{itemize}
'''



def parseArguments():
    parser = argparse.ArgumentParser(description='Produces Manhatten, QQ plot and supporting tex file  for output from some tools')
    parser.add_argument('--inp_asso',type=str,required=True, help="association files")
    parser.add_argument('--rsname',type=str,required=True, help="association files")
    parser.add_argument('--pheno',type=str,required=True,help="pheno names")
    parser.add_argument('--cov',type=str,help="pheno names", default="")
    parser.add_argument('--out', type=str,help="out of tex file",default="test.tex")
    parser.add_argument('--chro_header',type=str,required=True,help="chro header in inp files")
    parser.add_argument('--pos_header',type=str,required=True,help="pos header in inp files")
    parser.add_argument('--rs_header',type=str,required=True,help="rs header in inp files")
    parser.add_argument('--pval_header',type=str,required=True,help="pvalue header in inp files")
    parser.add_argument('--freq_header',type=str,required=False,help="pvalue header in inp files")
    parser.add_argument('--beta_header',type=str,required=True,help="beta header in inp files")
    parser.add_argument('--geno_plot',type=str,required=True,help="info program to be print in tex files")
    parser.add_argument('--locuszoom_plot',type=str,required=True,help="info program to be print in tex files")
    parser.add_argument('--annot_pdf',type=str,required=True,help="info program to be print in tex files")
    args = parser.parse_args()
    return args

args = parseArguments()
inp     = args.inp_asso
pheno= args.pheno
chro_ent=args.chro_header
pos_ent=args.pos_header
rs_ent=args.rs_header
beta_ent=args.beta_header
freq_ent=args.freq_header
pval_ent=args.pval_header
rsname=args.rsname
cov=args.cov

result = pd.read_csv(inp,delim_whitespace=True)
if freq_ent ==None :
   result[freq_ent]=-1
best_row = r"%s & %s & %s & %7.4f & %6.3E & %7.4f\\"
x=result.head(1)
rsvalue = best_row%(x[chro_ent][0],x[rs_ent][0].replace("_","-"),x[pos_ent][0],x[beta_ent][0],x[pval_ent][0], x[freq_ent][0])+EOL
tabres=r"""
\begin{table}[tb]
\begin{center}
\begin{tabular}{|l l l l r l|}\hline
Chr & SNP & Pos & Beta & P & Freq\\\hline
%(best)s\hline
\end{tabular}
\end{center}
\caption{gwas result for rs}
\label{tab:top:%(rsname)s}
\end{table}
\clearpage
"""

EOL = chr(10)
latex_tex=tabres
latex_tex+=r"""
\includepdf[pages=-]{"""+args.annot_pdf+"""}
"""

latex_tex+=r"""
\clearpage
\section{Locus Zoom plot}
\begin{figure}[ht]
\begin{center}
\includegraphics[width=15cm]{"""+args.locuszoom_plot.replace(".pdf","") + """}
\caption{%(rsname)s : locus zoom plot}
\label{fig:locuszoom}
\end{center}
\end{figure}
"""


latex_tex+=r"""
\clearpage
\section{Phenotype in function of genotypes}
\begin{figure}[ht]
\begin{center}
\includegraphics[width=15cm]{"""+args.geno_plot.replace(".pdf","")+"""}
\caption{%(rsname)s : Phenotype in function of genotypes}
\label{fig:genotypeplot}
\end{center}
\end{figure}
"""




latex_tex+="\end{document}\n"
hashd = { 'pheno':pheno.replace("_","-"), 'cov':cov.replace("_","-"), 'best':rsvalue, 'rsname':rsname, "date":check_output("date").strip()}
template=template.replace("*-",chr(92)).replace("##",chr(36)).replace("@.@",chr(10))
tex_file = template+latex_tex
tex_file=tex_file%hashd

g = open(args.out,"w")
g.write(tex_file)
g.close()



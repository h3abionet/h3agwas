#!/usr/bin/env python3
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import sys
import pandas as pd
import os
import argparse
import numpy as np
from  subprocess import CalledProcessError
import subprocess
import re

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

try:
   kpsewhich=check_output("which kpsewhich",shell=True)
   if kpsewhich:
      kpsewhich=check_output("kpsewhich datetime.sty",shell=True)
except CalledProcessError:
   kpsewhich=""

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
*-title{transformation pipeline}
*-date{%(date)s}
'''+dateheader+('''
*-author{H3Agwas : transformation of vcf in plink format}

*-newcommand{*-ourfig}[3]{*-begin{figure}[ht]*-begin{center}*-includegraphics[scale=0.6]{#3} *-end{center} *-caption{#2 [File is *-protect*-url{#3}]}  *-label{#1}*-end{figure}}
*-begin{document}
*-maketitle

*-section{Introduction}
This report gives a brief overview of transformation vcf file from imputation in plink format.
''')


def obtainstatfile(filel):
   readf=open(filel)
   linetest=readf.readline().split()
   readf.close()

   def splitfreq(splline):
     return(float(splline[5]))

   def cmpfreq(splline):
       return(float(spll[6])/float(spll[5]))

   if len(linetest)==6:
      funcfreq=splitfreq
   elif len(linetest)==7:
      funcfreq=cmpfreq
   else :
      print(linetest)
      print("error for computed frequencies " +str(len(linetest))+"\n")
      sys.exit(23)
      
   listfreq=[]
   listscore=[]
   readf=open(filel)
   for line in readf :
      spll=line.split()
      listfreq.append(funcfreq(spll))
      listscore.append(float(spll[4]))
   return (listfreq, listscore) 

def parseArguments():
    parser = argparse.ArgumentParser(description='extract annotation for specific position')
    parser.add_argument('--list_files',type=str,required=True, help="file contains chro and list of files with annotation", default=None)
    parser.add_argument('--min_score',type=float,required=True,help="file contain one position : rs chro pos")
    parser.add_argument('--out', type=str,help="output")
    args = parser.parse_args()
    return args
args = parseArguments()
TAB =chr(9)
lspl=args.list_files.split(',')
minscore=args.min_score

out=args.out

writestatfile=open(out+'.chro.stat', 'w')
writestatfile.write('File\tnbfreqnonnull\tnbscoreupper\tnbscorehighfreqnonnull\ttotal\n')

allstatfreq=[]
allstatscore=[]
statresume=[]
for filel in lspl :
   (listfreq, listscore)=obtainstatfile(filel)
   allstatfreq+=listfreq
   allstatscore+=listscore
   nbnonnullfreq=0
   nbscoregood=0
   nbnonnullfreqscoregood=0
   for cmt in range(len(listfreq)) :
     if listfreq[cmt]>0:
       nbnonnullfreq+=1  
     if listscore[cmt]>minscore:
       nbscoregood+=1
     if listfreq[cmt]>0 and listscore[cmt]>minscore :
       nbnonnullfreqscoregood+=1 
   #nbnonnullfreq=len([x for x in listfreq if x>0])
   #nbscoregood=len([x for x in listscore if x>minscore])
   nball=len(listfreq)
   statchro=[filel,nbnonnullfreq,nbscoregood,nbnonnullfreqscoregood,nball]
   writestatfile.write("\t".join([str(x) for x in statchro])+'\n')
   statresume.append(statchro)

nbnonnullfreq=0
nbscoregood=0
nbnonnullfreqscoregood=0
for cmt in range(len(allstatfreq)) :
   if allstatfreq[cmt]>0:
     nbnonnullfreq+=1
   if allstatscore[cmt]>minscore:
     nbscoregood+=1
   if allstatfreq[cmt]>0 and allstatscore[cmt]>minscore :
     nbnonnullfreqscoregood+=1
nball=len(listfreq)
statchro=["Total",nbnonnullfreq,nbscoregood,nbnonnullfreqscoregood,nball]
writestatfile.write("\t".join([str(x) for x in statchro])+'\n')
statresume.append(statchro)



def gettab(result, resume_row):
    best = ""
    for r in result:
        print(r)
        row = resume_row%("*-protect*-url{"+r[0]+"}", r[1], r[2], r[3], r[4])+EOL
        best = best + row
    return best


plt.hist(allstatscore,bins=10,range=(0,1),align="mid",rwidth=0.9,color="b",edgecolor="red")
plt.title("histogram of score")
plt.xlabel("Score")
plt.ylabel("Frequence")
plt.legend()
plt.show()
filescore=out.replace('.','_')+"_score.pdf"
plt.savefig(filescore)


plt.hist(allstatfreq,bins=10,range=(0,1),align="mid",rwidth=0.9,color="b",edgecolor="red")
plt.title("histogram of frequencies")
plt.xlabel("Frequencie")
plt.ylabel("Frequence")
plt.legend()
plt.show()
filefreq=out.replace('.','_')+"_freq.pdf"
plt.savefig(filefreq)

latex_tex = r"""
\clearpage
Resume by files format, SNPs number higher than the score higher than %(score)s and with frequencies more than 0
\begin{table}[tb]
\begin{center}
\begin{tabular}{|l|l|l|l|l|}\hline
File & $f > 0$ & $Sc > %(score)s$ & $Sc > %(score)s and f > 0$ & Total \\\hline
%(stats)s\hline
\end{tabular}
\end{center}
\caption{Table contains resume of score number (Sc) higher than upper limit and frequencies (f) not null by file }
\label{tab:top:resumefreqscore}
\end{table}

\noindent
The score distribution plot can be found in Figure \ref{fig:scoredistribution}.The frequencies distribution total
be found in Figure~\ref{fig:frequenciesdistribution}. 

\begin{figure}[ht]
\begin{center}
\includegraphics[width=15cm]{%(plotscore)s}
\caption{score distribution}
\label{fig:scoredistribution}
\end{center}
\end{figure}

\begin{figure}[ht]
\begin{center}
\includegraphics[width=15cm]{%(plotfreq)s}
\caption{Distribution of frequencies}
\label{fig:frequenciesdistribution}
\end{center}
\end{figure}

\clearpage
"""

resume_row = r"%s & %d & %d & %d & %d \\"
resumerow=gettab(statresume, resume_row)
hashd = {'stats':resumerow, 'score':args.min_score, 'plotfreq':filefreq,'plotscore':filescore, "date":check_output("date").strip()}

template=(template+ latex_tex)%hashd




template=template+r"""
*-end{document}
"""
template=template.replace("*-",chr(92)).replace("##",chr(36)).replace("@.@",chr(10))

outpdf=out+".tex"
g=open(outpdf,"w")
g.write(template)
g.close()

os.system("pdflatex %s >& /dev/null"%outpdf)
os.system("pdflatex %s"%out)


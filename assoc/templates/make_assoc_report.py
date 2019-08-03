#!/usr/bin/env python3

# (c) University of the Witwatersand, Johannesburg on behalf of the H3ABioinformatics Network Consortium
# 2016-2018
# Licensed under the MIT Licence
# See the "LICENSE" file for details

import glob
import sys
import os
from  subprocess import CalledProcessError
import subprocess
import re

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

if len(sys.argv)<=1:
   sys.argv = "make_assoc_report.py ${params.pheno} $texf".split()

pheno  = "*-protect*-url{%s}"%sys.argv[1]
out    = sys.argv[2]



src_tex  = sorted(glob.glob("*tex"))





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
*-title{Association Testing  %(base)s : %(pheno)s}
*-date{%(date)s}
'''+dateheader+'''
*-author{H3Agwas Association Testing Pipeline}

*-newcommand{*-ourfig}[3]{*-begin{figure}[ht]*-begin{center}*-includegraphics[scale=0.6]{#3} *-end{center} *-caption{#2 [File is *-protect*-url{#3}]}  *-label{#1}*-end{figure}}
*-begin{document}

*-maketitle

*-section{Introduction}

This report gives a brief overview of the run of the association testing pipeline.
*-begin{itemize}
*-item You were testing for the following phenotypes *-url{${these_phenos}}
*-item You were using the following covariates [*-url{${these_covariates}}]
*-end{itemize}
'''

EOL = chr(10)

for  x in src_tex:
    template = template + (r"*-input %s"%x) + EOL+"*-clearpage"+EOL

template = template + r"""

*-section{Technical details}

The analysis and report was produced by the h3aGWAS pipeline (*-url{http://github.com/h3abionet/h3agwas}) produced by the Pan-African Bioinformatics Network for H3Africa (*-url{http://www.h3abionet.org}).

The following tools were used:

*-begin{itemize}
*-item $nextflowversion [Di Tommaso et al, 2017]
*-item $wflowversion
*-item The command line below was called [NB: if the command line is long, the linebreak may break oddly after a hyphen or dash so take care.]
*-begin{lstlisting}
${workflow.commandLine}
*-end{lstlisting}
*-item The profile ${workflow.profile} was used%(dockerimages)s
*-item The full configuration can be found in the appendix.
*-end{itemize}


*-clearpage

*-appendix

*-section{Configuration}

The Nextflow configuration files are shown below.

%(config)s

*-end{document}
"""


pdict = {}
pdict["date"]=check_output("date").strip()
pdict["base"]="Draft"
pdict["pheno"]=pheno
config = """

$config

"""


pdict['config']=config

def getImages(images):
   images =images.replace("[","").replace("]","").replace(",","").split()
   result = "Table "+chr(92)+"ref{table:docker}"+chr(10)+chr(92)+"begin{table}"+chr(92)+"begin{tabular}{ll}"+chr(92)+"textbf{Nextflow process} &" + chr(92)+"textbf{Docker Image}"+chr(92)+chr(92)+chr(92)+"hline"+chr(10)
   for img in images:
      dets = img.split(":",1)
      if len(dets)==1:
         (proc,dimg)=("default",dets[0])
      else:
         (proc,dimg)==img.split(":",1)
      result = result +  \
                  proc.lstrip(chr(36)) + "&" + chr(92) + "url{%s}"%dimg+\
                  chr(92)+chr(92)
   result = result+chr(92)+"end{tabular}"+chr(92)+"caption{Docker Images Used}" +chr(92)+"label{table:docker}"+chr(92)+"end{table}"
   return result


images = "${workflow.container}"
if images=="[:]":
   pdict["dockerimages"] = ": locally installed binaries used"
else:
   images = getImages("${workflow.container}")
   pdict["dockerimages"] = ": the docker images used are found in "+images


template = template % pdict



template=template.replace("*-",chr(92)).replace("##",chr(36)).replace("@.@",chr(10))







g=open(out,"w")
g.write(template)
g.close()

os.system("pdflatex %s >& /dev/null"%out)
os.system("pdflatex %s"%out)

#!/usr/bin/env python3

# Scott Hazelhurst, 2016
# Creates a PDF report for QC
# 
# (c) University of the Witwatersand, Johannesburg on behalf of the H3ABioinformatics Network Consortium
# 2016-2018
# Licensed under the MIT Licence.
# See the "LICENSE" file for details
# 


from __future__ import print_function

import argparse
import sys
import os
import re
from  subprocess import CalledProcessError
import subprocess

def check_output(x,shell=False):
   ans=subprocess.check_output(x,shell=shell)
   ans=str(ans,'ascii') # encoding option only frmo python 3.6
   return ans


# Check we have pdflatex and up to date style file
#pdflatex=check_output("which pdflatex",shell=True)
#if len(pdflatex)<2:
#   check_output("echo no pdflatex > report.pdf",shell=True)
#  print("No PDFLATEX")
#   sys.exit(0)

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


# check if we are being called from the command line or as a template
# from Nextflow. If  as a template from Nextflow, then the nextflow process
# must define BASE
# NB: only for a Nextflow template -- you could also directly call from a Nextflow
# script if this script is on the path -- in which case the parameters should be passed
if len(sys.argv)==1:
   sys.argv = [sys.argv[0],"$orig","$base","$cbim","$cfam","$missingvhetpdf","$mafpdf","$dupf","$fsex","$indmisspdf","$snpmisspdf"]


def parseArguments():
    parser=argparse.ArgumentParser()
    parser.add_argument('orig', type=str),
    parser.add_argument('base', type=str),
    parser.add_argument('cbim', metavar='CBIM', type=str,help='clean bim'),
    parser.add_argument('cfam', metavar='CFAM', type=str,help='cleanfam'),
    parser.add_argument('missingvhetpdf', type=str),
    parser.add_argument('mafpdf', type=str),
    parser.add_argument('dupf', type=str),
    parser.add_argument('fsex', type=str),
    parser.add_argument('indmisspdf', type=str),
    parser.add_argument('snpmisspdf', type=str),        
    args = parser.parse_args()
    return args

args = parseArguments()
pdict = vars(args)


diff_miss_text_exists = """
Figure~*-ref{fig:diffP} plots the differences between
cases and controls, showing the SNP-wise p-value, unadjusted for multiple testing

*-ourfig{fig:diffP}{The plot shows differential missingness for each
  (log) level of significance, the number of SNPs with p-value of at
  least this significance. The flatter the better.}{$diffmisspdf}

For removal of SNPs, we compute the p-value adjusted for multiple
testing, by performing permutation testing (1000 rounds) using the
PLINK mperm maxT option.  SNPs are removed from the data set if their
adjusted (EMP2) differential missingness p-value is less than
${params.cut_diff_miss}. The SNPs that are removed can be found in the
file *-url{$diffmiss}
"""

diff_miss_text_not_exists = """
There are no SNPs with significant differential missingess in this data set. If your data has been very stringently QCd so
that there are zero or very few SNPs with *-emph{any} missingness, this is likely. However, in large data with a reasonable amount
of missing data  this is statistically unlikely and may indicate that something has gone wrong.
"""

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
*-title{Quality control report for *-protect*-url{%(base)s}}
*-date{%(date)s}
'''+dateheader+'''
*-author{H3Agwas QC Pipeline}

*-newcommand{*-ourfig}[3]{*-begin{figure}[ht]*-begin{center}*-includegraphics[scale=0.6]{#3} *-end{center} *-caption{#2 [File is *-protect*-url{#3}]}  *-label{#1}*-end{figure}}
*-begin{document}

*-maketitle

*-section{Introduction}

The input file for this analysis was *-url{%(orig)s.{bed,bim,fam}}. This data includes:
*-begin{itemize}
*-item %(numrsnps)s SNPs
*-item %(numrfam)s  participants
*-end{itemize}

*-noindent The input files and md5 sums were
*-begin{lstlisting}
%(inpmd5)s
*-end{lstlisting}

*-noindent
Note that some statistics are shown twice -- on the raw input data and on the final result, since these statistics are needed or different purposed.

*-subsection*{Approach}

The pipeline takes an incremental approach to QC, trading extra
computation time in order to achieve high quality while removing as
few data as possible. Rather than applying all cut-offs at once, we
incrementally apply cutoffs (for example, removing really badly
genotyped SNPs before checking for heterozygosity will result in fewer individuals
failing heterozygosity checks).


*-section{QC Phase 0}

*-label{sec:phasezero}

This phase only removes SNPs which are duplicated (based on SNP name). No other QC is done and so the output of this phase should really be considered as raw data.
*-begin{enumerate}
*-item There were %(numdups)s duplicate SNPs. The file with them (if any) is called *-protect*-url{%(dupf)s}. Note that duplicate SNPs are determined by the names of the SNPs. SNPs which appear at the same position are probably duplicates but may not be. You can control whether you want to detect these using the parameter *-url{remove_on_bp}. *-textbf{It is crucial to examine this file to avoid inadvertently removing SNPs. On some chips there are duplicate SNPs at a position -- you should select what you what want.}

*-item %(numfailedsex)s individuals had discordant sex information -- the full PLINK report can be found in *-url{%(fullsex)s} and an extract of the PLINK report showing only the failed reports can be found can be found in *-url{%(fsex)s}, and a more detailed analysis can be found in Section~*-ref{sec:batch}.
*-end{enumerate}

Figure *-ref{fig:snpmiss} shows the spread of missingness per SNP across the sample, whereas Figure *-ref{fig:indmiss} shows the spread of missingness per individual across the sample. Note that this shows missingness before any filtering or cleaning up of the data.

*-ourfig{fig:snpmiss}{SNP missingness: For each level of missingness specified on the ##x## axis, the corresponding ##y##-value shows the proportion of SNPs which have missingness *-emph{less} than this.}{{%(snpmisspdf)s}}

*-ourfig{fig:indmiss}{Missingness per indvididual: For each level of missingness specified on the ##x## axis, the corresponding ##y##-value shows the proportion of individuals which have missingness *-emph{less} than this.}{{%(indmisspdf)s}}


*-input $initmaftex

*-input $inithwetex

*-clearpage
*-section{QC Phase 1}

The details of the final filtering can be found in the Nextflow script. Note that the exact ordering or removal will affect the final results. However, we take a conservative approach.


*-begin{enumerate}
*-item Only autosomal SNPs are included
*-item SNPs and individuals that have been very poorly genotyped (missingness exceeding 10 per cent) are removed.
*-item Individuals with missingness greater than ${params.cut_mind} are removed (by filtering out very badly genotyped individuals or SNPs) in the previous steps we *-emph{may} save a few individuals in this step;
*-item SNP with missingness greater than ${params.cut_geno} are removed ;
*-item minor allele frequency less than ${params.cut_maf} (and greater than 1-${params.cut_maf});
*-item HWE p-value less than ${params.cut_hwe}
*-end{enumerate}

*-noindent

*-input $qc1



*-input $batch_tex

*-clearpage


*-section{Heterozygosity check}

Levels of heterozygosity were examined in the data filtered in the
previous step. Figure~*-ref{fig:missvhet} shows plots of
heterozygosity versus individual missingness (i.e., the number of SNPs
missing per individual).  Levels of heterozygosity should be between
the ranges given in the confguration file -- anything higher may indicate that there is sample
contamination, lower may indicate inbreeding. However, each set of
data must be treated on its own merits and the analyst must apply
their mind the problem. Missingness should be low.

*-ourfig{fig:missvhet}{Missingness versus heterozygosity: the lines show the mean heterozygosity plus/minus standard deviations (the brighter the colours, the greater the density). If there is zero missingness, only heterozygosity is shown in a violin  plot.}{{%(missingvhetpdf)s}}

Individuals out of range heterozygosity were removed. Any indviduals with heterozygosity:

*-begin{itemize}
*-item less than ${params.cut_het_low} are removed. This may indicate inbreeding.
*-item greater than ${params.cut_het_high} are removed. This may indicate sample contamination.
*-end{itemize}
 
*-noindent
Overall %(numhetrem)s individuals were removed. These individuals, if any, can be found in the file *-url{$misshetremf}.




*-clearpage
*-section{Minor Allele Frequency Spread}

Figure~*-ref{fig:maf} shows the cumulative distribution of minor
allele frequency in the data *-textbf{after} quality control (the figures shown in Section *-ref{sec:phasezero} show the MAF before QC. The MAF cut-off should be chosen high enough that one is sure that the variants seen are real (so this would depend on the size of the sample and the quality of the genotyping and whether some of the data is imputed). In this analysis the cut off was ${params.cut_maf}. Again, note that the *-emph{minor} allele is determined with respect to the frequency spectrum in this data -- `minor' is not synonym for alternate or non-reference allele, or the allele that has minor frequency in some other data set. Under this definition the MAF is always ##*-leq 0.5##.


*-ourfig{fig:maf}{Minor allele frequency distribution}{{%(mafpdf)s}}

*-clearpage

*-section{Differences between cases and controls}

We do not expect there to be large, observable macro-scale differences between cases and controls. Great caution needs to be taken in this case. If the samples are from heterogeneous groups, where the case-control status differs between groups, then there may well statistically significant differences between the cases and controls and a batch analysis should be undertaken. 

We compute for each SNP the missingness in the cases, and the
missingness in the controls, and the corresponding p-value describing
the difference in missingness. 

We expect very few SNPs to have highly significant differences. Where
many SNPs with very highly significant p-values are found, great care
should be taken. 

%(diff_miss_text)s


Figure~*-ref{fig:pca} shows a principal component analysis of the
data, identifying the cases and controls. Should the cases and
controls cluster differ signficantly, something is likely wrong.
Moreover should there be any significant clusters or outliers, association
testing should take into account stratification. Statistical testing could also
be done. Figure~*-ref{fig:pcaeigen} shows for each principal component what the corresponding eigenvalue is. The shape of the curve may indicate if there is any structure in the data. Informally, the ``broken stick" model sees two processes generating PCs -- population structure and random fluctuation. The eigenvalue is a function of both -- the former has a major impact but drops relatively rapidly; the latter has less effect but drops more slowly. If (!) the model is correct then you will see an inflection point in the graph, where the random effects become the leading factor in the change.  More formal analysis may be desirable, e.g., using the Tracy-Widom statistic or Velicer's MAP test.

*-ourfig{fig:pca}{Principal Component Analysis of Cases Versus Controls}{{$pcapdf}}

*-ourfig{fig:pcaeigen}{Eigenvalues for each principal component: the shape of the curve gives some indication of how many PCs are important}{{$eigenvalpdf}}

*-clearpage
*-section{Hardy-Weinberg Equilibrium}

Deviation for Hardy-Weinberg Equilibrium (HWE) may indicate sample contamination. However, this need not apply to cases, nor in a situation where there is admixture. For each SNP, we compute  the probability of the null hypothesis (that  the deviation from HWE is by chance alone).  Figure~*-ref{fig:hwe} shows a plot of the corresponding p-value versus the frequency of occurrence.

*-ourfig{fig:hwe}{The plot shows for each level of significance, the number of SNPs with H
WE p-value}{{$hwepdf}}






*-section{Final results}

*-noindent
The quality control procedures as described below were applied to the data. 
The final, cleaned result contains:
*-begin{itemize}
*-item %(numcsnps)s SNPs
*-item %(numcfam)s  participants
*-end{itemize}

*-noindent
The final output files are 
*-begin{itemize}
*-item *-url{$cbed}, 
*-item *-url{$cbim}, and 
*-item *-url{$cfam}.
*-end{itemize}

The ouput files' md5 sums are shown below

*-begin{verbatim}
%(outmd5)s
*-end{verbatim}

*-pagebreak[4]

*-section{Technical details}

The analysis and report was produced by the h3aGWAS pipeline (*-url{http://github.com/h3abionet/h3agwas}) produced by the Pan-African Bioinformatics Network for H3Africa (*-url{http://www.h3abionet.org}).

The following tools were used:

*-begin{itemize}
*-item %(plinkversion)s  [Chang et al 2015]
*-item $nextflowversion [Di Tommaso et al, 2017]
*-item $wflowversion
*-item The command line below was called [NB: if the command line is long, the linebreak may break oddly after a hyphen or dash so take care.]
*-begin{lstlisting}
${workflow.commandLine}
*-end{lstlisting}
*-item The profile ${workflow.profile} was used%(dockerimages)s
*-item The full configuration can be found in the appendix.
*-end{itemize}

*-pagebreak[4]
*-section{References}

*-begin{itemize}
*-item Chang  CC, Chow CC, Tellier LC, Vattikuti S, Purcell SM, and Lee JJ. (2015). Second-generation PLINK: rising to the challenge of larger and richer datasets. *-emph{GigaScience}, 4(1), 1-16. *-url{http://doi.org/10.1186/s13742-015-0047-8}
*-item Di Tommaso P, Chatzou M, Prieto Barja P, Palumbo P, Notredame C. Nextflow enables reproducible computational workflows. *-emph{Nature Biotechnology}, 35, 316-319, 2017. Nextflow can be downloaded from *-url{https://www.nextflow.io/}
*-end{itemize}

*-clearpage



*-appendix
*-section{Configuration}

{*-footnotesize


%(configuration)s


}
*-end{document}'''

# gymnastics to get backslash and dollar
template=template.replace("*-",chr(92)).replace("##",chr(36))

def countLines(fn):
    count=0
    with open(fn) as f:
        for  line in f:
            count=count+1
    return count

def readLines(fn):
    resp=""
    with open(fn) as f:
        for  line in f:
            resp=resp+line
    return resp

def getImages(images):
   images =images.replace("[","").replace("]","").replace(",","").split()
   result = "Table "+chr(92)+"ref{table:docker}"+chr(10)+chr(92)+"begin{table}"+chr(92)+"begin{tabular}{ll}"+chr(92)+"textbf{Nextflow process} &" + chr(92)+"textbf{Docker Image}"+chr(92)+chr(92)+chr(92)+"hline"+chr(10)
   for img in images:
      dets = img.split(":",1)
      if len(dets)==1:
         (proc,dimg)=("default",dets[0])
      else:
         (proc,dimg)=dets
      result = result +  \
                  proc + "&" + chr(92) + "url{%s}"%dimg+\
                  chr(92)+chr(92)
   result = result+chr(92)+"end{tabular}"+chr(92)+"caption{Docker Images Used}" +chr(92)+"label{table:docker}"+chr(92)+"end{table}"
   return result


s=os.stat("$diffmisspdf")

diff_miss_text = diff_miss_text_exists if s.st_size>0 else diff_miss_text_not_exists
pdict["diff_miss_text"] = diff_miss_text
 
f=open(args.orig)
pdict['numrsnps'] = f.readline().split()[0]
pdict['numrfam']  = f.readline().split()[0]
f.close()

f=open("$ilog")
pdict['plinkversion']=f.readline()
f.close()


pdict['numhetrem'] =  countLines("$misshetremf")
pdict['numcsnps'] =  countLines(args.cbim)
pdict['numcfam']  =  countLines(args.cfam)
pdict['numdups']  =  countLines(args.dupf)
pdict['numdiffmiss'] = countLines("$diffmiss")
pdict['irem']       = "$irem"

pdict['inpmd5']= readLines("$inpmd5")
pdict['outmd5']= readLines("$outmd5")


pdict['configuration']="""$config_text""".replace("*-",chr(92)).replace("##",chr(36)).replace("@.@",chr(10))

if "${workflow.container}"=="[:]":
   pdict["dockerimages"] = ": locally installed binaries used"
else:
   images = getImages("${workflow.container}")
   pdict["dockerimages"] = ": the docker images used are found in "+images


pdict["dockerimages"]=pdict["dockerimages"].replace(chr(36),"")

# changed to support Python 3.5 and before
pdict["date"]=check_output("date").strip()


num_fs = countLines(args.fsex)
if num_fs == 1:
    head=open(args.fsex).readline()
    if "No sex" in head: num_fs=0
pdict['numfailedsex']=num_fs
pdict['fullsex']=re.sub("badsex","sexcheck",args.fsex)

out=open("%s.tex"%args.base,"w")
text = template%pdict
text=text.replace("*-",chr(92)).replace("##",chr(36))
out.write (text)
out.close()

os.system("pdflatex %s >& /dev/null"%args.base)
os.system("pdflatex %s"%args.base)
  

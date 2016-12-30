#!/usr/bin/env python

# Scott Hazelhurst, 2016
# Creates a PDF report for QC

# Tested under both Python 2.7 and 3.5.2
# 
# Scott Hazelhurst on behalf of the H3ABioinformatics Network Consortium
# December 2016
# (c) Released under GPL v.2


from __future__ import print_function

import argparse
import sys
import os
import re
from  subprocess import check_output, CalledProcessError


# Check we have pdflatex and up to date style file
pdflatex=check_output("which pdflatex",shell=True)
if len(pdflatex)<2:
   check_output("echo no pdflatex > report.pdf",shell=True)
   print("No PDFLATEX")
   sys.exit(0)

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



template='''
*-documentclass[11pt]{article}

*-usepackage[paper=a4paper,left=2cm,right=2cm,top=2cm,bottom=2cm]{geometry}
*-usepackage{graphicx}
*-usepackage{url}
*-title{Quality control report for %(base)s}
*-date{%(date)s}
'''+dateheader+'''
*-author{H3Agwas QC Pipeline}

*-newcommand{*-ourfig}[3]{*-begin{figure}[ht]*-begin{center}*-includegraphics[scale=0.6]{#3} *-end{center} *-caption{#2 [File is #3]}  *-label{#1}*-end{figure}}
*-begin{document}

*-maketitle

*-section{Introduction}

The input file for this analysis was *-url{%(base)s.{bed,bim,fam}}. This data includes:
*-begin{itemize}
*-item %(numrsnps)s SNPs
*-item %(numrfam)s  participants
*-end{itemize}

*-noindent The input files and md5 sums were
*-begin{verbatim}
%(inpmd5)s
*-end{verbatim}

*-subsection*{Approach}

The pipeline takes an incremental approach to QC, trading extra
computation time in order to achieve high quality while removing as
few data as possible. Rather than applying all cut-offs at once, we
incrementally apply cutoffs (for example, removing really badly
genotyped SNPs before checking for heterozygosity wil result in fewer individuals
failing heterozygosity checks).


*-section{Initial QC}

*-begin{enumerate}
*-item There were %(numdups)s duplicate SNPs. The file with them (if any) is called {*-em %(dupf)s}.
*-item %(numfailedsex)s individuals had discordant sex information -- further information can be found in {*-em %(fsex)s}.
*-end{enumerate}


*-section{Heterozygosity check}

Levels of heterozygosity were examined. Figure~*-ref{fig:missvhet} shows plots of heterozygosity versus
individual missingness (i.e., the number of SNPs missing per
individual).  Levels of heterozygosity should be between the ranges
given -- anything higher may indicate that there is sample
contamination, lower may indicate inbreeding. However, you need to
apply your mind to the data. Missingness should be low.

*-ourfig{fig:missvhet}{Missingness versus heterozygosity}{*-detokenize{%(missingvhetpdf)s}}

Individuals out of range heterozygosity were removed: first
SNPs with a genotyping failure rate of over 10 per cent were removed,
and then heterozysgosity checked. Any indviduals with heterozygosity:

*-begin{itemize}
*-item less than ${params.cut_het_low} are removed. This may indicate inbreeding.
*-item greater than ${params.cut_het_high} are removed. This may indicate sample contamination.
*-end{itemize}
 
*-noindent
Overall %(numhetrem)s individuals were removed. These individuals, if any, can be found in the file *-url{$misshetremf}.


Figure *-ref{fig:snpmiss} shows the spread of missingness per SNP across the sample, whereas *-ref{fig:indmiss} show the spread of missingness per individual across the sample these should be compared.

*-ourfig{fig:snpmiss}{SNP missingness}{*-detokenize{%(snpmisspdf)s}}

*-ourfig{fig:indmiss}{Missingness per indvididual}{*-detokenize{%(indmisspdf)s}}






*-pagebreak[3]
*-section{Minor Allele Frequency Spread}

Figure~*-ref{fig:maf} shows the cumulative distribution of minor
allele frequency in the data. The MAF cut-off should be chosen high enough that you are sure that the variants you are seeing are real (so this would depend on the size of the sample). You have chosen a cut off of ${params.cut_maf}.


*-ourfig{fig:maf}{Minor allele frequency distribution}{*-detokenize{%(mafpdf)s}}


*-section{Relatedness}

Using PLINK, relatedness is computed using IBD and ##{*-widehat{*-pi}}## as a proxy. All pairs of individuals with a ##{*-widehat{*-pi} *-geq ${pi_hat} }## are examined -- that individual with the greater missingness is removed. The ##*-widehat{*-pi}## of ${pi_hat} is a parameter of the pipeline. 

%(numrels)s individuals were removed because of relatedness. The list of all individuals can be found in the *-url{$relf} file. 

*-section{Differences between cases and controls}

We do not expect there to be large, observable macro-scale differences between cases and controls. Great caution needs to be taken in this case. 

We compute for each SNP the missingness in the cases, and the
missingness in the controls, and the corresponding p-value describing
the difference in missingness. 

We expect very few SNPs to have highly significant differences. Where
many SNPs with very highly significant p-values are found, great care
should be taken. Figure~*-ref{fig:diffP} plots the differences between
cases and controls, showing the SNP-wise p-value, unadjusted for multiple testing

*-ourfig{fig:diffP}{The plot shows for each (log) level of significance, the number of SNPs with that p-value}{*-detokenize{$diffmisspdf}}

For removal of SNPs, we compute the p-value adjusted for multiple testing, by performing permutation testing (1000 rounds) using the PLINK mperm maxT option.
SNPs are removed from the data set if their adjusted (EMP2) differential missingness p-value is less than ${params.cut_diff_miss}. The SNPs that are removed can be
found in the file *-url{$diffmiss}



Figure~*-ref{fig:pca} shows a principal component analysis of the
data, identifying the cases and controls. Should the cases and
controls cluster differ signficantly, something is likely wrong.
Moreover should there be any significant clusters or outliers, association
testing should take into account stratification. Statistical testing could also
be done.

*-ourfig{fig:pca}{Principal Component Analysis of Cases Versus Controls}{*-detokenize{$pcapdf}}

*-section{Hardy-Weinberg Equilibrium}

Deviation for Hardy-Weinberg Equilibrium (HWE) may indicate sample contamination. However, this need not apply to cases, nor in a situation where there is admixture. For each SNP, we compute  the probability of the null hypothesis (that  the deviation from HWE is by chance alone).  Figure~*-ref{fig:hwe} shows a plot of the corresponding p-value versus the frequency of occurrence.

*-ourfig{fig:hwe}{The plot shows for each level of significance, the number of SNPs with H
WE p-value}{*-detokenize{$hwepdf}}

The description of how SNPs were filtered based on HWE are discussed in Section~*-ref{sec:final}


*-section{Final filtering}

*-label{sec:final}
The details of the final filtering can be found in the Nextflow script. Note that the exact ordering or removal will affect the final results. However, we take a conservative approach.


*-begin{enumerate}
*-item SNPs that failed differential missingness,  and individuals that have been very poorly genotyped (missingness exceeding 20 per cent) are removed.
*-item Then, SNPs that have been very poorly genotyped (missingness exceeding 20 per cent) are removed.
*-item Finally we select only autosomal SNPs and filter out SNPs  with
*-begin{itemize}
*-item minor allele frequence less than ${params.cut_maf};
*-item individual missingness greater than ${params.cut_mind};
*-item SNP missingness greater than ${params.cut_geno}; and 
*-item HWE p-value less than ${params.cut_hwe}
*-end{itemize}
*-end{enumerate}

The individuals removed in this phase, if any, can be found in the file *-url{%(irem)s}.

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
*-item R version %(rversion)s [R Core Team, 2016]
*-item $nextflowversion [Di Tommaso et al]
*-item $wflowversion
*-item The command line *-verb:${workflow.commandLine}: was called
*-item The profile ${workflow.profile} was used%(dockerimages)s
*-item The full configuration can be found in the appendix.
*-end{itemize}

*-pagebreak[4]
*-section{References}

*-begin{itemize}
*-item Chang, C. C., Chow, C. C., Tellier, L. C., Vattikuti, S., Purcell, S. M., and Lee, J. J. (2015). Second-generation PLINK: rising to the challenge of larger and richer datasets. *-emph{GigaScience}, 4(1), 1-16. *-url{http://doi.org/10.1186/s13742-015-0047-8}
*-item R Core Team (2016). *-emph{R: A language and environment for statistical
  computing}. R Foundation for Statistical Computing, Vienna, Austria.
  *-url{https://www.R-project.org/}
*- Paolo Di Tommaso, Maria Chatzou, Pablo Prieto Baraja, Cedric Notredame. A novel tool for highly scalable computational pipelines. *-url{http://dx.doi.org/10.6084/m9.figshare.1254958}. Nextflow can be downloaded from *-url{https://www.nextflow.io/}
*-end{itemize}

*-clearpage

*-appendix
*-section{nextflow.config}

{*-footnotesize

*-begin{verbatim}
%(configuration)s
*-end{verbatim}

}
*-end{document}'''

# gymnastics to get backslash and dollar
template=template.replace("*-",unichr(92)).replace("##",unichr(36))

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
   result = "Table "+unichr(92)+"ref{table:docker}"+unichr(10)+unichr(92)+"begin{table}"+unichr(92)+"begin{tabular}{ll}"+unichr(92)+"textbf{Nextflow process} &" + unichr(92)+"textbf{Docker Image}"+unichr(92)+unichr(92)+unichr(92)+"hline"+unichr(10)
   for img in images:
      (proc,dimg)=img.split(":")
      result = result +  \
                  proc + "&" + unichr(92) + "url{%s}"%dimg+\
                  unichr(92)+unichr(92)
   result = result+unichr(92)+"end{tabular}"+unichr(92)+"caption{Docker Images Used}" +unichr(92)+"label{table:docker}"+unichr(92)+"end{table}"
   return result
 
f=open(args.orig)
pdict['numrsnps'] = f.readline().split()[0]
pdict['numrfam']  = f.readline().split()[0]
f.close()

f=open("$ilog")
pdict['plinkversion']=f.readline()
f.close()

f=open("rversion")
pdict['rversion']=f.readline()
f.close()

pdict['numhetrem'] =  countLines("$misshetremf")
pdict['numcsnps'] =  countLines(args.cbim)
pdict['numcfam']  =  countLines(args.cfam)
pdict['numdups']  =  countLines(args.dupf)
pdict['numdiffmiss'] = countLines("$diffmiss")
pdict['numrels']       = countLines("$relf")
pdict['irem']       = "$irem"

pdict['inpmd5']= readLines("$inpmd5")
pdict['outmd5']= readLines("$outmd5")


f=open("$nextflowconfig")
conf=""
for line in f: conf=conf+line
f.close()
pdict['configuration']=conf

if "${workflow.container}"=="[:]":
   pdict["dockerimages"] = ": locally installed binaries used"
else:
   images = getImages("${workflow.container}")
   pdict["dockerimages"] = ": the docker images used are found in "+images

pdict["dockerimages"]=pdict["dockerimages"].replace(unichr(36),"")
pdict["date"]=check_output("date")

num_fs = countLines(args.fsex)
if num_fs == 1:
    head=open(args.fsex).readline()
    if "No sex" in head: num_fs=0
pdict['numfailedsex']=num_fs

    
out=open("%s.tex"%args.base,"w")
out.write (template%pdict)
out.close()

os.system("pdflatex %s >& /dev/null"%args.base)
os.system("pdflatex %s"%args.base)
  

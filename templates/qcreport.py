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

*-usepackage{a4wide}
*-usepackage{graphicx}
*-usepackage{url}

*-title{Quality control report for %(base)s}

*-author{H3Agwas QC Pipeline}

*-newcommand{*-ourfig}[3]{*-begin{figure}[ht]*-begin{center}*-includegraphics[scale=0.6]{#3} *-end{center} *-caption{#2 [File is #3]}  *-label{#1}*-end{figure}}
*-begin{document}

*-maketitle

*-section{Introduction}

The input file for this analysis was *-emph{%(base)s}. This data includes:
*-begin{itemize}
*-item %(numrsnps)s SNPs
*-item %(numrfam)s  participants
*-end{itemize}

*-noindent
The final, cleaned result contains:
*-begin{itemize}
*-item %(numcsnps)s SNPs
*-item %(numcfam)s  participants
*-end{itemize}


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








*-noindent
Figure *-ref{fig:snpmiss} shows the spread of missingness per SNP across the sample, whereas *-ref{fig:indmiss} show the spread of missingness per individual across the sample these should be compared.

*-ourfig{fig:snpmiss}{SNP missingness}{*-detokenize{%(snpmisspdf)s}}

*-ourfig{fig:indmiss}{Missingness per indvididual}{*-detokenize{%(indmisspdf)s}}







*-section{Minor Allele Frequency Spread}

Figure~*-ref{fig:maf} shows the cumulative distribution of minor
allele frequency in the data. The MAF cut-off should be chosen high enough that you are sure that the variants you are seeing are real (so this would depend on the size of the sample). You have chosen a cut off of ${params.cut_maf}.


*-ourfig{fig:maf}{Minor allele frequency distribution}{*-detokenize{%(mafpdf)s}}


*-section{Differences between cases and controls}

We do not expect there to be large, observable macro-scale differences between cases and controls. Great caution needs to be taken in this case. 

We compute for each SNP the missingness in the cases, and the
missingness in the controls, and the corresponding p-value
describing the difference in missingness. We expect very few SNPs to
have highly significant differences. Where many SNPs with very highly significant p-values are 
found, great care should be taken. Figure~*-ref{fig:diffP} plots the differences between cases and controls.

*-ourfig{fig:diffP}{The plot shows for each (log) level of significance, the number of SNPs with that p-value}{*-detokenize{$diffmisspdf}}




*-end{document}'''

template=template.replace("*-",unichr(92))

def countLines(fn):
    count=0
    with open(fn) as f:
        for  line in f:
            count=count+1
    return count


f=open(args.orig)
pdict['numrsnps'] = f.readline().rstrip()
pdict['numrfam']  = f.readline().rstrip()
f.close()



pdict['numhetrem'] =  countLines("$misshetremf")
pdict['numcsnps'] =  countLines(args.cbim)
pdict['numcfam']  =  countLines(args.cfam)
pdict['numdups']  =  countLines(args.dupf)

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
  

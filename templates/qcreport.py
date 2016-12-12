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
    sys.argv = [sys.argv[0],"$base","$rbim","$rfam","$cbim","$cfam","$missingvhetpdf","$mafpdf","$duplog","$dupf","$fsex","$indmiss","$snpmiss"]


def parseArguments():
    parser=argparse.ArgumentParser()
    parser.add_argument('base', metavar='BASE', type=str)
    parser.add_argument('rbim', metavar='RBIM', type=str,help='raw bim'),
    parser.add_argument('rfam', metavar='RFAM', type=str,help='raw fam'),
    parser.add_argument('cbim', metavar='CBIM', type=str,help='clean bim'),
    parser.add_argument('cfam', metavar='CFAM', type=str,help='cleanfam'),
    parser.add_argument('missingvhetpdf', type=str),
    parser.add_argument('mafpdf', type=str),
    parser.add_argument('duplog', type=str),
    parser.add_argument('dupf', type=str),
    parser.add_argument('fsex', type=str),
    parser.add_argument('indmiss', type=str),
    parser.add_argument('snpmiss', type=str),        
    args = parser.parse_args()
    return args

args = parseArguments()
pdict = vars(args)

template='''
*-documentclass[11pt]{article}

*-usepackage{a4wide}
*-usepackage{graphicx}

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

*-noindent
Individuals with out of range heterozygosity were removed. For this phase, SNPs with a genotyping failure rate of over 10 per cent were removed, and then heterozysgosity checked. Any indviduals with heterozygosity:

*-begin{itemize}
*-item less than ${params.cut_het_low} are removed. This may indicate inbreeding.
*-item greater than ${params.cut_het_high} are removed. This may indicate sample contamination.
*-end{itemize}
 
*-noindent
Overall %(numhetrem)s individuals were removed. These individuals, if any, can be found in the file *-emph{%(hetrem)s}.








*-noindent
Figure *-ref{fig:snpmiss} shows the spread of missingness per SNP across the sample, whereas *-ref{fig:indmiss} show the spread of missingness per individual across the sample these should be compared.

*-ourfig{fig:snpmiss}{SNP missingness}{*-detokenize{%(snpmiss)s}}

*-ourfig{fig:indmiss}{Missingness per indvididual}{*-detokenize{%(indmiss)s}}







*-subsection{Minor Allele Frequency Spread}

Figure~*-ref{fig:maf} shows the cumulative distribution of minor
allele frequency in the data.


*-ourfig{fig:maf}{Minor allele frequency distribution}{*-detokenize{%(mafpdf)s}}

*-end{document}'''

template=template.replace("*-",unichr(92))

def countLines(fn):
    count=0
    with open(fn) as f:
        for  line in f:
            count=count+1
    return count




pdict['numrsnps'] =  countLines(args.rbim)
pdict['numrfam']  =  countLines(args.rfam)
pdict['numdups']  =  countLines(args.dupf)
pdict['numcsnps'] =  countLines(args.cbim)
pdict['numcfam']  =  countLines(args.cfam)


num_fs = countLines(args.fsex)
if num_fs == 1:
    head=open(args.fsex).readline()
    if "No sex" in head: num_fs=0

pdict['numfailedsex']=num_fs
    
out=open("%s.tex"%args.base,"w")
out.write (template%pdict)
out.close()
os.system("pdflatex %s"%args.base)
os.system("pdflatex %s"%args.base)
  

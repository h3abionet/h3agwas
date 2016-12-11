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
    sys.argv = [sys.argv[0],"$base","$rbim","$rfam","$missingvhetpdf","$mafpdf","$duplog","$dupf","$fsex"]


def parseArguments():
    parser=argparse.ArgumentParser()
    parser.add_argument('base', metavar='BASE', type=str)
    parser.add_argument('rbim', metavar='RBIM', type=str,help='raw bim'),
    parser.add_argument('rfam', metavar='RFAM', type=str,help='raw fam'),
    parser.add_argument('missingvhetpdf', type=str),
    parser.add_argument('mafpdf', type=str),
    parser.add_argument('duplog', type=str),
    parser.add_argument('dupf', type=str),
    parser.add_argument('fsex', type=str),    
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

The input file for this analysis was {base}. This data includes:
*-begin{itemize}
*-item %(numrsnps)s SNPs
*-item %(numrfam)s  participants
*-end{itemize}

*-section{Quality control}

*-begin{enumerate}
*-item There were %(numdups)s duplicate SNPs. The file with them (if any) is called {*-em %(dupf)s}.
*-item %(numfailedsex)s individuals had discordant sex information -- further information can be found in {*-em %(fsex)s}.
*-end{enumerate}

*-subsection{Missingness versus heterozygosity}

Figure~*-ref{fig:missvhet} shows plots of heterozygosity versus
missingness.  Levels of heterozygosity should be between the ranges
given -- anything higher may indicate that there is sample
contamination, lower may indicate inbreeding. However, you need to
apply your mind to the data. Missingness should be good.

*-ourfig{fig:missvhet}{Missingness versus heterozygosity}{*-detokenize{%(missingvhetpdf)s}}


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




pdict['numrsnps']=countLines(args.rbim)
pdict['numrfam']=countLines(args.rfam)
pdict['numdups'] =  countLines(args.dupf)

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

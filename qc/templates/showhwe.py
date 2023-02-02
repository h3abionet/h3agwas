#!/usr/bin/env python3

import pandas as pd
import numpy  as np
import sys
from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt
import matplotlib
EOL=chr(10)


if len(sys.argv)<=1:
    sys.argv="showmaf.py $hwe  $base".split()



template = """
*-paragraph*{Hardy Weinberg Statistics}: 
Figure~*-ref{fig:hweinit} shows the cumulative distribution of Hardy-Weinberg ##p##-value for  the SNPs in the raw data. This can be used to assess the cost of excluding SNPs with a particular ##p##-cutoff. We expect the curve the fit tightly to the main diagonal, except for a very small ##p## values (and this deviation may not be observable on a linear plot).

*-begin{figure}[htb]*-centering
*-includegraphics[width=10cm]{%s}
*-caption{HWE distribution. For an HWE-value shown on the ##x##-axis, the corresponding ##y##-value shows the proportion of SNPs with HWE p-value *-emph{at least this frequency}; that is, it shows the proportion of SNPs which will *-emph{be removed} if the HWE filter of this ##x##-value is chosen. File is *-protect*-url{%s}.}
*-label{fig:hweinit}
*-end{figure}
"""

# We don't draw a Manhatten plot at the moment -- would be interesting to see if there are regions
# that are problematic -- problem is that the .hwe file doesn't have base position so we'd need to change
# the pipeline but it would be good to do
def drawManhatten(pheno, result,outf):
    fig, (ax1, ax2) =  plt.subplots(2, 1, sharey=True)
    ax=ax1
    delta=0
    colours= ['crimson','blue','green']
    xtick_pos=[]
    xtick_label = []
    chroms = result.groupby("CHR")
    for chrom_num, chrom_res in chroms:
        this_chrom = result['CHR']==chrom_num
        result.loc[this_chrom,'BP']+=delta
        old_delta=delta
        delta = delta + int(chrom_res.tail(1)['BP'])
        xtick_pos.append((delta+old_delta)/2)
        xtick_label.append(str(chrom_num))
        under_thresh = result['P']<0.005
        ax.scatter(result.loc[this_chrom & under_thresh, 'BP'],\
                   -np.log10(result.loc[this_chrom  & under_thresh,'P']),c=colours[chrom_num%3])
        if chrom_num == 9:
           ax.set_xticklabels(xtick_label)
           ax.set_xticks(xtick_pos)
           xtick_pos=[]
           xtick_label=[]
           ax=ax2
           delta=0
    ax.set_xticklabels(xtick_label)
    ax.set_xticks(xtick_pos)
    plt.savefig(outf)

def drawQQ(result,outf):
    plt.figure()
    fig, ax = plt.subplots()
    sort_p = -np.log10(result['P'].sort_values())
    n=len(sort_p)
    expected = -np.log10(np.linspace(1/n,1,n))
    ax.set_xlabel("HWE expected (-log p)",fontsize=14)
    ax.set_ylabel("HWE observed (-log p)",fontsize=14)
    plt.plot(expected,sort_p)
    plt.plot(expected,expected)
    plt.savefig(outf)


qq_template = """

The QQ plot for the HWE scores can be found in Figure *-ref{fig:qq}. The region of deviation from the line of expected versus observed ##p##-values will be more observable here. Note that if there are very small observed ##p##-values in relation to expected values, the expected curve may be very flat --- pay attention to the x and y axis coordinates.  Since we are plotting on a negative log-scale, note that regions of low probabiliy of deviation from HWE (##p##-value close to 1) are at the left, and regions of high probability (low ##p##-value) are at the right. The tail of the plot where deviation from the diagonal occurs is likely to be a good cut-off to use for QC.  

However, care needs to be taken not exclude SNPs. We are using HWE ##p##-value as a proxy for something having gone wrong with the sample or genotyping, and this is a little crude. In a study with participants from different population groups in a recently admixed group, deviation from HWE is expected and does not indicate problems with QC. Moreover, in a disease study, it is likely that those individuals that are affected, those SNPs that are associated with the condition under study will not in be in HWE. Care needs to be taken -- it is easier to handle in a pure case/control study. In a population cross-section study with different conditions being considered, it might be advisable to re-run the QC pipeline for HWE for each study.  The current version of the pipeline does not support his more complex analysis, though we plan to extend. 

*-begin{figure}[ht]
*-begin{center}
*-includegraphics[width=10cm]{%(qqfile)s}
*-end{center}
*-caption{QQ plot for Hardy-Weinberg scores. The graphic can be found in the file *-protect*-url{%(qqfile)s}}
*-label{fig:qq}
*-end{figure}

"""



def rename(x):
    if x.left < 0:
       return "0"
    else:
       return str(x)  




def getPic(frm,test,pdfout):
   fig = plt.figure(figsize=(17,14))
   fig,ax = plt.subplots()
   matplotlib.rcParams['ytick.labelsize']=13
   matplotlib.rcParams['xtick.labelsize']=13
   hwe = frm[frm["TEST"]==test]["P"]
   big = min(hwe.mean()+2*hwe.std(),hwe.nlargest(4).iloc[3])
   if big > max(0.95*len(hwe),100):
       hwe = hwe[hwe<big]
   hwe = np.sort(hwe)
   n = np.arange(1,len(hwe)+1) / np.float64(len(hwe))
   ax.step(hwe,n)
   ax.set_xlabel("HWE p score",fontsize=14)
   ax.set_ylabel("Proportion of SNPs with HWE p-value or less",fontsize=14)
   ax.set_title("Cumulative prop. of SNPs with HWE or less",fontsize=16)
   fig.tight_layout()
   plt.savefig(pdfout)


f = open(sys.argv[1])
header=f.readline()
test=False
for i in range(5):
    line=f.readline()
    if "ALL(QT)" in line:
        test="ALL(QT)"
    elif "ALL(NP)" in line:
        test="ALL(NP)"
    elif "ALL" in line:
        test="ALL"
if not test:
    print((EOL*5)+"The Hardy-Weinberg file is malformed, can't find ALL test <%s>"%sys.argv[1])
    sys.exit(12)

frm = pd.read_csv(sys.argv[1],delim_whitespace=True)
frm = frm[frm["TEST"]==test]
base = sys.argv[2]
pdfout = "%s.pdf"%base
texout = "%s.tex"%base
qqpdf  = "%s-qq.pdf"%base
getPic(frm,test,pdfout)
drawQQ(frm, qqpdf)

g=open(texout,"w")
g.write((template%(pdfout,pdfout)).replace("*-",chr(92)).replace("##",chr(36)))
g.write((qq_template%({'qqfile': qqpdf})).replace("*-",chr(92)).replace("##",chr(36)))
g.close()

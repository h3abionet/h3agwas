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
Figure~*-ref{fig:hweinit} shows the cumulative distribution of Hardy-Weinberg ##p##-value for  the SNPs in the raw data. This can be used to assess the cost of excluding SNPs with a particular ##p##-cutoff.

*-begin{figure}[htb]*-centering
*-includegraphics[width=10cm]{%s}
*-caption{HWE distribution. For an HWE-value shown on the ##x##-axis, the corresponding ##y##-value shows the proportion of SNPs with HWE p-value *-emph{at least this frequency}; that is, it shows the proportion of SNPs which will *-emph{be removed} if the HWE filter of this ##x##-value is chosen. File is *-protect*-url{%s}.}
*-label{fig:hweinit}
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
   hwe = np.sort(hwe[hwe<big])
   n = np.arange(1,len(hwe)+1) / np.float(len(hwe))
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

getPic(frm,test,pdfout)

g=open(texout,"w")
g.write((template%(pdfout,pdfout)).replace("*-",chr(92)).replace("##",chr(36)))
g.close()

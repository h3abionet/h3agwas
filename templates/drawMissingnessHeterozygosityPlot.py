#!//usr/bin/env python3

# (c) University of the Witwatersand, Johannesburg on behalf of the H3ABioinformatics Network Consortium
# 2016-2018
# Licensed under the MIT Licence
# See the "LICENSE" file for details

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator
import pandas as pd
import numpy as np
import argparse
import sys
import math
import bisect
from scipy.stats import gaussian_kde

colour_choices=["black","magenta","darkcyan","red","blue","orange"]

def parseArguments():
   if len(sys.argv)<=1:
      sys.argv=\
      "missHetPlot.py $cohortImiss $cohortHet $output".split()
   parser=argparse.ArgumentParser()
   parser.add_argument('imiss', type=str, metavar='imiss'),
   parser.add_argument('het', type=str, metavar='het'),
   parser.add_argument('output', type=str, metavar='output'),
   args = parser.parse_args()
   return args



args = parseArguments()
imiss  = pd.read_csv(args.imiss,delim_whitespace=True,usecols=["FID","IID","F_MISS"],index_col=["FID","IID"])
het    = pd.read_csv(args.het,delim_whitespace=True,index_col=["FID","IID"])



fig, ax = plt.subplots(figsize=(8,6))
font = {'family' : 'normal','weight' : 'bold','size'   : 14}
matplotlib.rc('font', **font)
matplotlib.rcParams['xtick.labelsize']=13
matplotlib.rcParams['ytick.labelsize']=13
matplotlib.rcParams['axes.titlesize']=16

het["meanHet"] = (het["N(NM)"] - het["O(HOM)"])/het["N(NM)"]
het["meanHet"] = np.where(np.isnan(het["meanHet"]), 0,het["meanHet"])
het_min = het["meanHet"].mean()-3*het["meanHet"].std()
het_max = het["meanHet"].mean()+3*het["meanHet"].std()
plt_min = max(0,het_min-3*het["meanHet"].std(),het["meanHet"].min())
plt_max = min(1,het_max+3*het["meanHet"].std(),het["meanHet"].max())
y=het["meanHet"]
#NOTE: IF F_MISS IS ZERO THEN WE ONLY PLOT MEAN HETEROZYGOSITY
if not  (imiss["F_MISS"].min() == imiss["F_MISS"].max() == 0):
   imiss["CALL_RATE"]= 1-imiss["F_MISS"]
   imiss["logF_MISS"] = np.where(imiss["F_MISS"]==0,0.00000001,np.log10(imiss["F_MISS"]))
   comb = het.merge(imiss,left_index=True,right_index=True)
   x  = imiss["logF_MISS"]
   xy = np.vstack([x,y])
   dens_colsx = gaussian_kde(xy)(xy)

   ax.scatter(x,y,c=dens_colsx,s=1)
   plt.title("Samples plotted by missingness versus heterozygosity")
   plt.xlabel("Missingness genotypes",fontsize=14)
   plt.ylabel("Heterozygosity rate",fontsize=14)

   ax.set_xticks([-3,-2,-1])
   ax.set_xticklabels([0.001,0.01,0.1])
   yticks=[0,0.05,0.10,0.15,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.4,0.45,0.5]
   start = bisect.bisect(yticks,plt_min)
   if yticks[start]>=plt_min : start = start-1
   end   = bisect.bisect(yticks,plt_max)+1
   yticks = yticks[start:end]
   #plt.ylim(yticks[0]-0.025,yticks[-1]+0.025)
   ax.set_yticks(yticks)
   ax.axhline(y=het_min,color="red",alpha=0.4)
   ax.axhline(y=het_max,color="blue",alpha=0.4)

else:
   ax.set_xticks([])
   ax.violinplot(y)
   plt.xlabel("Proportion of SNPs")
   plt.ylabel("Heterozygosity")
   plt.title("Violin plot showing distribution of heterozygosity")



plt.tight_layout()
plt.savefig(args.output,type="pdf")

#!/usr/bin/env python3


import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

inp = sys.argv[1]
out_man = sys.argv[2]
out_qq  = sys.argv[3]
result = pd.read_csv(inp,delim_whitespace=True)
chroms = result.groupby("chr")


fig, (ax1, ax2) =  plt.subplots(2, 1, sharey=True)

ax=ax1
delta=0
colours= ['crimson','blue','green']
xtick_pos=[]
xtick_label = []
for chrom_num, chrom_res in chroms:
    this_chrom = result['chr']==chrom_num
    result.loc[this_chrom,'ps']+=delta
    old_delta=delta
    delta = delta + int(chrom_res.tail(1)['ps'])
    xtick_pos.append((delta+old_delta)/2)
    xtick_label.append(str(chrom_num))
    under_thresh = result['p_wald']<0.005
    ax.scatter(result.loc[this_chrom & under_thresh, 'ps'],\
               -np.log10(result.loc[this_chrom  & under_thresh,'p_wald']),c=colours[chrom_num%3])
    if chrom_num == 9:
       ax.set_xticklabels(xtick_label)
       ax.set_xticks(xtick_pos)
       xtick_pos=[]
       xtick_label=[]
       ax=ax2
       delta=0
ax.set_xticklabels(xtick_label)
ax.set_xticks(xtick_pos)

plt.savefig(out_man)

plt.figure()
sort_p = -np.log10(result['p_wald'].sort_values())
n=len(sort_p)
expected = -np.log10(np.linspace(1/n,1,n))
print(n,len(expected))
plt.plot(expected,sort_p)
plt.plot(expected,expected)
plt.savefig(out_qq)


#!/usr/bin/env python3
#Load SNP frequency file and generate histogram

import pandas as   pd
import numpy  as np
import sys
from matplotlib import use
use('Agg')
import argparse
import matplotlib
import matplotlib.pyplot as plt
import sys

def parseArguments():
    if len(sys.argv)<=1:
        sys.argv="snpmissPlot.py $input $label $output".split()
    parser=argparse.ArgumentParser()
    parser.add_argument('input', type=str, metavar='input'),
    parser.add_argument('label', type=str, metavar='label'),
    parser.add_argument('output', type=str, metavar='output'),
    args = parser.parse_args()
    return args

args = parseArguments()

data = pd.read_csv(args.input,delim_whitespace=True)

fig = plt.figure(figsize=(17,14))
fig,ax = plt.subplots()
matplotlib.rcParams['ytick.labelsize']=13
matplotlib.rcParams['xtick.labelsize']=13
miss = data["F_MISS"]
big = min(miss.mean()+2*miss.std(),miss.nlargest(4).iloc[3])
interesting = miss[miss<big]
if len(interesting)>0.9 * len(miss):
    miss = interesting
miss = np.sort(miss)
n = np.arange(1,len(miss)+1) / np.float(len(miss))
ax.step(miss,n)
ax.set_xlabel("Missingness",fontsize=14)
ax.set_ylabel("Proportion of %s"%args.label,fontsize=14)
ax.set_title("Cumulative prop. of  %s with given missingness"%args.label,fontsize=16)
fig.tight_layout()
plt.savefig(args.output)

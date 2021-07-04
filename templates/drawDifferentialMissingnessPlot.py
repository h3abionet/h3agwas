#!/usr/bin/env python3

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

def parseArguments():
   if len(sys.argv)<=1:
      sys.argv=\
      "drawDifferentialMissingnessPlot.py $input $output".split() 
   parser=argparse.ArgumentParser()
   parser.add_argument('input', type=str, metavar='imiss'),
   parser.add_argument('output', type=str, metavar='output'),
   args = parser.parse_args()
   return args


args = parseArguments()
frq = pd.read_csv(args.input,delim_whitespace=True)
if len(frq) >= 1:
  frq["logP"] = np.log10(frq["P"])
  fig, ax = plt.subplots(figsize=(9,8))
  font = {'family' : 'normal','weight' : 'bold','size'   : 13}
  matplotlib.rc('font', **font)
  matplotlib.rcParams['xtick.labelsize']=15
  matplotlib.rcParams['ytick.labelsize']=15
  miss=np.sort(frq["logP"])[1:]
  n = np.arange(1,len(miss)+1) / np.float(len(miss))
  ax.step(miss,n)
  ax.set_xlabel("logP differential missingness",fontsize=14)
  ax.set_ylabel("Fraction of SNPs",fontsize=14)
  ax.set_title("Cumulative proportion of SNPs with given differential missingness.")
  fig.tight_layout()
  plt.savefig(args.output)
else:
    g=open(args.output,"w")
    g.close()

#!//usr/bin/env python3


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator
import pandas as pd
import numpy as np
import argparse
import sys

colour_choices=["black","magenta","darkcyan","red","blue","orange"]

def parseArguments():
   if len(sys.argv)<=1:
      sys.argv=\
      "drawPCA.py $base $cc $cc_fname $col $eigvals $eigvecs $output".split()
   parser=argparse.ArgumentParser()
   parser.add_argument('input', type=str, metavar='input'),
   parser.add_argument('cc', type=str, metavar='label'),
   parser.add_argument('cc_fname', type=str, metavar='label'),
   parser.add_argument('column', type=str, metavar='label'),
   parser.add_argument('eigvals', type=str, metavar='label'),
   parser.add_argument('eigvecs', type=str, metavar='output'),
   parser.add_argument('output', type=str, metavar='output'),
   args = parser.parse_args()
   return args



def getColours():
    if len(args.cc_fname)==0:
        return [], "black"
    phe = pd.read_csv(args.cc,delim_whitespace=True)
    def our_colour(x):
        return colour_choices[x[args.column]]
    the_colours = phe.apply(our_colour,axis=1)
    options = ["missing","control","case"]
    # There may be no missing items
    return list(enumerate(options))[phe[args.column].min():], the_colours


def getEigens():
    evals = np.loadtxt(args.eigvals)
    pc1 = 100*evals[0]/evals[:10].sum()
    pc2 = 100*evals[1]/evals[:10].sum()
    evecs  = pd.read_csv(args.eigvecs, delim_whitespace=True, header=None)
    return pc1, pc2, evecs

col_names=['FID','IID']+list(map(lambda x: "PC%d"%x,range(1,21)))


def draw(pc1,pc2,evecs,labels,the_colours):
   fig = plt.figure(figsize=(12,12))
   fig, ax = plt.subplots()
   locator = MaxNLocator(nbins=5) 
   ax.xaxis.set_major_locator(locator)
   ax.scatter(evecs[2],evecs[3],s=1,c=the_colours)
   ax.legend(scatterpoints=1)
   recs=[]
   classes=[]
   for (i,label) in labels:
      recs.append(mpatches.Rectangle((0,0),1,1,fc=colour_choices[i]))
      classes.append(label)
   plt.legend(recs,classes,loc=4)
   plt.xlabel("PC1 (variation %4.1f %%)"%pc1)
   plt.ylabel("PC2 (variation %4.1f %%)"%pc2)
   plt.tight_layout()
   plt.savefig(args.output,type="pdf")



args = parseArguments()
labels, the_colours =  getColours()
evecs =  getEigens()
draw(*evecs,labels,the_colours)

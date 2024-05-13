#! /usr/bin/env python3


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator
import pandas as pd
import numpy as np
import argparse
import sys

EOL=chr(10)
gap = EOL*3

colour_choices=["black","magenta","darkcyan","red","blue","orange","aqua","beige","chartreuse","darkblue","gold","indigo","ivory","olive","sienna","wheat","salmon","orangered","silver","tan","grey","lightblue","violet","yellow","turquoise", "yellowgreen","khaki","goldenrod","aquamarine","azure","brown","crimson","fuchsia"]

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
    if len(args.cc_fname)==0 or args.cc_fname=="0":
        return [], "_None_"
    phe = pd.read_csv(args.cc,delim_whitespace=True)
    all_labels = phe[args.column].unique()
    all_labels.sort()
    colours = colour_choices
    while len(all_labels) > len(colours):
       colours=colours+colour_choices
    the_colour_choices = dict(zip(all_labels,colours[:len(all_labels)]))
    def our_colour(x):
        try:
           result = the_colour_choices[x[args.column]]
        except:
           sys.exit(gap+"There's a problem with the phenotype file <%s>, column <%s>, ID <%s>%s"%(args.cc,args.column,x,gap))
        return result
    phe['colour'] = phe.apply(our_colour,axis=1)
    return list(enumerate(all_labels)), phe


def getEigens():
    evals = np.loadtxt(args.eigvals)
    fig, ax = plt.subplots(figsize=(10,8))
    ax.plot(range(len(evals)),evals)
    plt.xlabel("Principal component",fontsize=15)
    plt.ylabel("Eigenvalue",fontsize=15)
    plt.xticks(range(len(evals)+1), list(map(lambda x : x+1, range(len(evals))) )+[""])
    plt.tight_layout()
    plt.savefig("eigenvalue.pdf")
    pc1 = 100*evals[0]/evals[:10].sum()
    pc2 = 100*evals[1]/evals[:10].sum()
    pcs = list(map(lambda x: "PC%d"%(x+1),range(len(evals))))
    evecs  = pd.read_csv(args.eigvecs, delim_whitespace=True, header=None,names=["FID","IID"]+pcs)
    return pc1, pc2, evecs

col_names=['FID','IID']+list(map(lambda x: "PC%d"%x,range(1,21)))


def draw(pc1,pc2,evecs,labels,phe):
   if type(phe)==str and phe == "_None_":
       evecs["colour"] = "black"
   else:
       evecs=evecs.merge(phe,how='inner',on=["FID","IID"])   
   fig, ax = plt.subplots(figsize=(10,8))
   font = {'family' : 'normal','weight' : 'bold','size'   : 14}
   matplotlib.rc('font', **font)
   matplotlib.rcParams['xtick.labelsize']=13
   matplotlib.rcParams['ytick.labelsize']=13
   locator = MaxNLocator(nbins=5) 
   ax.xaxis.set_major_locator(locator)
   ax.scatter(evecs['PC1'],evecs['PC2'],s=1,c=evecs['colour'])
   ax.legend(scatterpoints=1)
   recs=[]
   classes=[]
   for (i,label) in labels:
      recs.append(mpatches.Rectangle((0,0),1,1,fc=colour_choices[i]))
      classes.append(label)
   plt.legend(recs,classes,loc=4)
   plt.xlabel("PC1 (variation %4.1f %%)"%pc1,fontsize=15)
   plt.ylabel("PC2 (variation %4.1f %%)"%pc2,fontsize=15)
   plt.tight_layout()
   plt.savefig(args.output)
   fig, ax = plt.subplots(figsize=(10,8))


def outputTeX():
     g=open("B040-pca.tex","w")
     g.write(("""
*-section{Principal Component Analysis of Participants}

Figure *-ref{fig:pca} shows a PCA of the participants. This should be examined for possible structure.

*-begin{figure}[htb]
*-centering
*-includegraphics[width=12cm]{%s}
*-caption{PCA of participants}
*-label{fig:pca}
*-end{figure}

     """%(args.output)).replace("*-",chr(92)).replace("##",chr(36)).replace("@.@",chr(10)))
     g.close()

args = parseArguments()
labels, phe =  getColours()
evecs =  getEigens()
draw(*evecs,labels,phe)
outputTeX()

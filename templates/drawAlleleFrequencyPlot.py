#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy  as np
import sys
import matplotlib
from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt


def main():
    EOL=chr(10)
    args=parseArguments()
    frm = pd.read_csv(args.input,delim_whitespace=True)
    getPic(frm,args)


def parseArguments():
    if len(sys.argv)<=1:
        sys.argv="drawAlleleFrequencyPlot.py $input $output".split()
    parser=argparse.ArgumentParser()
    parser.add_argument('input', type=str, metavar='input'),
    parser.add_argument('output', type=str, metavar='output'),
    args = parser.parse_args()
    return args

def getPic(frm,args):
    mafs = np.sort(frm['MAF'])
    n = np.arange(1,len(mafs)+1) / np.float(len(mafs))
    fig,ax = plt.subplots()
    ax.step(mafs,n)
    matplotlib.rcParams['xtick.labelsize']=13
    matplotlib.rcParams['ytick.labelsize']=13
    ax.set_xlabel("Minor allele frequency",fontsize=14)
    ax.set_ylabel("Proportion of SNPS",fontsize=14)
    plt.title('Cumulative MAF-spectrum on QC-ed data',fontsize=16)
    plt.savefig(args.output)


if __name__ == '__main__':
    main()
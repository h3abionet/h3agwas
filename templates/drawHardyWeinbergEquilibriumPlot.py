#!/usr/bin/env python3
#Load HWE P-value file and generate frequency_distribution

import pandas as   pd
import numpy  as np
import sys
from matplotlib import use
use('Agg')
import argparse
import matplotlib.pyplot as plt
import sys

def main():
    args = parseArguments()
    fig,ax = plt.subplots(figsize=(8,6))
    data = pd.read_csv(args.input,delim_whitespace=True)
    vals = np.sort(np.log10(data["P"]))[1:]
    n = np.arange(1,len(vals)+1) / np.float(len(vals))
    ax.step(vals,n)
    ax.set_xlabel("logP (HWE)")
    ax.set_ylabel("Proportion of SNPs")
    fig.tight_layout()
    plt.ylim(-0.01,1.01)
    plt.savefig(args.output)

def parseArguments():
    if len(sys.argv)<=1:
        sys.argv="drawHardyWeinbergEquilibriumPlot.py $input $output".split()
    parser=argparse.ArgumentParser()
    parser.add_argument('input', type=str, metavar='input'),
    parser.add_argument('output', type=str, metavar='output'),
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()
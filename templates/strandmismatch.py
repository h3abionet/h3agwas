#!/usr/bin/env python3

import argparse
import sys
import re



def parseArguments():
    parser=argparse.ArgumentParser()
    parser.add_argument('strandreport', type=str, metavar='strandreport'),
    parser.add_argument('flips', type=str, metavar='flips',help="flips"),
    args = parser.parse_args()
    return args


# If called with no parameters, we assume using Nexflow's template mechanism and the parameters
# are substituted in with fixed names
if len(sys.argv)<=1:
    sys.argv = ["strandmismatch.py","$strandreport","$flips"]


# we avoid backslashes
TAB=chr(9)
EOL=chr(10)

def getHeadings(f):
    line=f.readline()
    while line[0]=="#": 
        line=f.readline()
    data=line.split()
    print(data)
    col_f1   = data.index("Forward_Allele1")
    col_f2   = data.index("Forward_Allele2")
    col_topA = data.index("Top_AlleleA")
    col_topB = data.index("Top_AlleleB")
    col_snp  = data.index("SNP_Name")
    return [col_snp,col_topA, col_topB, col_f1, col_f2]


def produceFlips(f,col_snp,col_topA, col_topB, col_f1, col_f2):
    g = open(args.flips,"w")
    for line in f:
        try:
            data = line.rstrip().split()
            if data[col_topA] != data[col_f1]:
                g.write("{}{}".format(data[col_snp],EOL))
        except IndexError:
            pass
    g.close()

args = parseArguments()
if args.strandreport == "empty.txt":
    g=open(args.flips,"w")
    g.close()
    sys.exit(0)
f = open(args.strandreport)
cols = getHeadings(f)
produceFlips(f,*cols)

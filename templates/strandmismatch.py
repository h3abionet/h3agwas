#!/usr/bin/env python3

import argparse
import sys
import re
import pandas as pd
from operator import xor

def parseArguments():
    parser=argparse.ArgumentParser()
    parser.add_argument('strandreport', type=str, metavar='strandreport'),
    parser.add_argument('manifest', type=str, metavar=''),
    parser.add_argument('output_align', type=str, metavar='output_align'),
    parser.add_argument('flips', type=str, metavar='flips',help="flips"),
    args = parser.parse_args()
    return args


# If called with no parameters, we assume using Nexflow's template mechanism and the parameters
# are substituted in with fixed names
if len(sys.argv)<=1:
    sys.argv = ["strandmismatch.py","$strandreport","$manifest","$output_align","$flips"]


# we avoid backslashes
TAB=chr(9)
EOL=chr(10)

args = parseArguments()

g=open(args.flips,"w")
if args.output_align == "topbottom":
    g.close()
    sys.exit(0)


if args.output_align == "dbsnp":
    strand = pd.read_csv(args.strandreport,delim_whitespace=True,comment="#",\
                     dtype={"Chr":str},\
                     usecols=["SNP_Name","Chr","Coord","Top_AlleleA","Forward_Allele1"])
    snps_to_flip = strand[strand["Top_AlleleA"]!=strand["Forward_Allele1"]]["SNP_Name"]
elif args.output_align == "ref":
    mf = pd.read_csv(args.manifest,skiprows=7,delimiter=",",\
                     dtype={"Chr":str},\
                     usecols=["Name","IlmnStrand","RefStrand","Chr","MapInfo"])
    snps_to_flip = mf[(mf["IlmnStrand"]=="BOT")^(mf["RefStrand"]=="-")]["Name"]
else:
    g.close()
    sys.exit(args.output_align+" is not a supported output format")

g.write(EOL.join(snps_to_flip.values))
g.close()


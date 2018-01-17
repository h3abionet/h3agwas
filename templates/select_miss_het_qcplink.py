#!/usr/bin/env python

# This script selects those who fail heterozygosity constraints -- previous versions also looked at missingness too
# hence the name


from pandas import read_csv
import sys

TAB = chr(9)
if len(sys.argv)<=1:
  sys.argv = ["select_miss_het.qcplink.py","$het",float("${params.cut_het_low}"),float("${params.cut_het_high}"),"$outfname"]

hetf = sys.argv[1]
cut_het_high=float(sys.argv[3])
cut_het_low=float(sys.argv[2]);
outfname = sys.argv[4]


het      = read_csv(hetf,delim_whitespace=True)
mean_het = (het["N(NM)"]-het["O(HOM)"])/het["N(NM)"]
failed   = het[(mean_het<cut_het_low) | (mean_het>cut_het_high)]

failed.to_csv(outfname,header=False,sep=TAB)

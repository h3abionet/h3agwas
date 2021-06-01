#!/usr/bin/env python
import datatable as dt
import pandas as pd
import sys, os.path

# create file to write out lgen data
gsgtFile = sys.argv[1]
threads = sys.argv[2]

outFile = gsgtFile.replace(".csv.gz",".lgen")
lgenFile = open(os.path.basename(outFile), "w")

# load genotype reports and select necessary columns
gsgtFrame = dt.fread(gsgtFile, nthreads=threads, skip_to_line=11, na_strings=["-"])
lgenFrameDatatable = gsgtFrame[:, {"FID": gsgtFrame[1], "IID": gsgtFrame[1], "RSID": gsgtFrame[0], "Allele1": gsgtFrame[2], "Allele2": gsgtFrame[3]}]

# convert datatable frame to pandas frame
lgenFramePandas = lgenFrameDatatable.to_pandas()
lgenFramePandas.to_csv(lgenFile, sep=" ", na_rep="0", header=False, index=False, index_label=False, chunksize=1000)


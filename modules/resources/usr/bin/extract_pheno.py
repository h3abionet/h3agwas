#!/usr/bin/env python3
# (c) University of the Witwatersand, Johannesburg on behalf of the H3ABioinformatics Network Consortium
# 2016-2018
# Licensed under the Creative Commons Attribution 4.0 International Licence. 
# See the "LICENSE" file for details

import sys
import pandas as pd
import numpy as np

EOL=chr(10)

def check_missing(x):
    if x in ["NA","na","null","-9"]  or ((type(x)==type("x")) and (len(x)<1)) or x==-9 or x==r"\N" or x=="-":
        return "NA"
    else:
        return x

def errorMessage10(phe):
    print("""

    A problem has been detected in file <%s> column <%s>.

    There is some invalid data. I regret I can't tell you which row.


    Please check -- the data should be numeric only.


    If there is missing data, please use   NA



    """%(sys.argv[1],phe))


dataf = pd.read_csv(sys.argv[1],delim_whitespace=True)
columns = dataf.columns

pheno_labels_0 = sys.argv[2].split(",")


pheno_labels = ["FID","IID"]
if "FID" not in columns or "IID" not in columns:
    sys.exit(EOL*3+"It is mandatory for FID IID to be columns of the file %s"%sys.argv[1]+EOL*3)

for lab in pheno_labels_0:
    det = lab.split("/")
    phe = det[0]
    if phe not in columns:
        sys.exit((EOL*3+"<%s> given as phenotype, but is not a column of the file <%s>"+EOL+EOL)%(phe,sys.argv[1]))
    if len(det)>1:
        fn = eval(det[1])
        try:
            dataf[det[0]]=fn(dataf[phe])
        except:
            errorMessage10(phe)
            sys.exit(10)
    pheno_labels.append(phe)
    dataf[phe]=dataf[phe].apply(check_missing)


dataf[pheno_labels].to_csv(sys.argv[3],na_rep="-9",sep=chr(9),index=False,header=True)

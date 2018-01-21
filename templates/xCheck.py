#!/usr/bin/env python3


import argparse
import sys
import re
import os
import pandas as pd
import numpy as np
if len(sys.argv)<=1:
    sys.argv=["xCheck.py","$x","$out"]
    missingness=eval("$missingness")
else:
    missingness = [0.01,0.03,0.05]


def getResForM(base,m):
    out = "%s-%s"%(base,m)
    os.system("plink --bfile %s --mind %s --check-sex --out %s"%(base,m,out))
    sf = pd.read_csv("%s.sexcheck"%out,delim_whitespace=True,index_col=[0,1],usecols=["FID","IID","STATUS","F"])
    sf[m] = np.where(sf['STATUS']=='OK',"OK", np.where((sf['F']>=0.34) & (sf['F']<=0.66),"S","H"))
    return sf


base   = sys.argv[1]
outfn  = sys.argv[2]

result = getResForM(base,1)
for m in missingness:
    try:
        curr   = getResForM(base,m)
    except IOError:
        result[m]=np.nan
        continue
    result = result.join(curr[m],how='outer')
result=result.drop(['STATUS','F'],axis=1)
pd.to_pickle(result,outfn)



#!/usr/bin/env python3

import sys
import pandas as pd

EOL = chr(10)

if len(sys.argv) <= 1:
    sys.argv=["select_diffmiss_qcplink.py","$missing", "$probcol", "$cuff_diff_miss","$failed"]

mfr     = pd.read_csv(sys.argv[1],delim_whitespace=True)
probcol = sys.argv[2]
cut_diff_miss = float(sys.argv[3])

wanted = mfr[mfr[probcol]<cut_diff_miss]["SNP"].values
out = open(sys.argv[4],"w")
out.write(EOL.join(map(str,wanted)))
out.close()

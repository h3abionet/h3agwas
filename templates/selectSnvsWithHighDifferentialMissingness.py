#!/usr/bin/env python3

import sys
import pandas as pd

EOL = chr(10)

if len(sys.argv) <= 1:
    sys.argv=[
        "selectSamplesWithHighDifferentialMissingness.py",
        "$input",
        "$probcol",
        "${params.differentialMissingness.cut}",
        "$output"]
probcol = sys.argv[2]
cut_diff_miss = float(sys.argv[3])

mfr = pd.read_csv("$input",delim_whitespace=True)
wanted = mfr[mfr["EMP2"]<${params.differentialMissingness.cut}]["SNP"].values
out = open($output,"w")
out.write(EOL.join(map(str,wanted)))
out.close()

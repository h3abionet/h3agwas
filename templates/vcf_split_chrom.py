#!/usr/bin/env python3

import sys
import re
import gzip
import os
if len(sys.argv)==1:
    sys.argv=["vcf_split_chrom.py","$inputf"]



fname =sys.argv[1]
f = gzip.open(fname)


fname_f = os.path.basename(fname)
m = re.search("(.*).vcf.gz",fname_f)
base = m.group(1)

header = ""

line=""
while not re.search("#CHROM.*POS",line):
    line = f.readline()
    header=header+line

g=False

regexp = re.compile("^("+unichr(92)+"w+).*")
oldchrom="******"
for line in f:
    m = regexp.match(line)
    chrom = m.group(1)
    if g and chrom != oldchrom:
        g.close()
        g=False
    if not g:
        g = gzip.open("{}-{}.vcf.gz".format(base,chrom),"w")
        g.write(header)
    g.write(line)
    oldchrom=chrom
f.close()

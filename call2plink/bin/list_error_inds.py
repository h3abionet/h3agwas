#!/usr/bin/env python3

import sys
import re

TAB=chr(9)
EOL=chr(10)


mask_type = sys.argv[1]
mask_file = sys.argv[2]
id_pat    = sys.argv[3]
fam       = sys.argv[4]
out       = sys.argv[5]


if mask_type == "fid-iid":
    def transform(line):
        m= re.search(id_pat, line)
        if m:
            fid, iid=m.group(1),m.group(2)
        else:
            sys.exit("Can't extract FID, IID with <%s> from <%s>"%(id_pat,line))
        return (fid,iid)
elif mask_type == "sample-label":
    def transform(line):
        data = line.strip()
        return (data[0],data[0])
else:
    sys.exit("Mask_type <%s> is incorrect: can only be fid-iid or sample-label")


dels = []
for line in open(mask_file):
    print(len(line))
    data = line.split()
    print(data[0].strip())
    if mask_type=="sample-label":
        dels.append((data[0].strip(),data[0].strip()))
    else:
        dels.append((data[0].rstrip(),data[1].rstrip()))

f = open(fam)
g=open(out,"w")
for line in f:
    data = line.split()
    fid = data[0].strip()
    iid = data[1].strip()
    real_id = transform(fid)
    if real_id in dels:
        g.write("%s %s%s"%(fid,iid,EOL))
g.close()

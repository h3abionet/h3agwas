#!/usr/bin/env python3

import sys
import re

TAB=chr(9)
EOL=chr(10)

fname = sys.argv[1]
gname = sys.argv[2]


f = open(fname)
g = open(gname,"w")

IDS=set()
for line in f:
    data = line.split()
    the_id = data[0]+":"+data[1]
    if the_id in IDS:
        m=re.search("(.*DUP)(\d+)",data[1])
        if m:
            data[1] = m.group(1)+str(int(m.group(2))+1)
        else:
            data[1] = data[1]+"DUP{}".format(0)
    else:
        IDS.add(the_id)
    g.write(TAB.join(data)+EOL)
g.close()


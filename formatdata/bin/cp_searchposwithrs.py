#!/usr/bin/env python3

import sys
from itertools import chain
FileRsI=open(sys.argv[1])
FileAllRs=sys.stdin
WriteRs=open(sys.argv[2],'w')
PosRs=int(sys.argv[3])


ListeRs=set(list(chain(*[x.replace('\n','').split(';') for x in FileRsI if len(x)>1])))
for l in FileAllRs :
    if l[0]=="#" :
       continue
    ls=l.split('\t')
    if len(ls)>PosRs and ls[PosRs] in ListeRs :
      WriteRs.write(l)

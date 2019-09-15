#!/usr/bin/env python3

# extract from bim file duplicate bim

import sys

readBim=open(sys.argv[1])
readRs=open(sys.argv[2])
writenewbim=open(sys.argv[3],'w')
ListRs=set([x.replace('\n','') for x in readRs])
for line in readBim :
    tmp=line.split()
    if tmp[1] in ListRs or tmp[1]=='.':
        tmp[1]=tmp[0]+':'+tmp[3]
    writenewbim.write("\t".join(tmp)+'\n')

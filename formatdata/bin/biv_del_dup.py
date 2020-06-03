#!/usr/bin/env python3
import sys

filei=open(sys.argv[1])
fileout=open(sys.argv[2], 'w')
listrs=set([])
for x in filei :
   splx=x.split()
   if splx[0] not in listrs :
      fileout.write(x)
      listrs.add(splx[0])

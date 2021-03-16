#!/usr/bin/env python3
import sys
readf=open(sys.argv[1])
perc=float(sys.argv[2])
nbind=0

for line in readf :
  nbind+=1

readf.close()
controlind=int(nbind*(1-perc))
caseind=nbind -controlind
print(caseind," ",controlind)

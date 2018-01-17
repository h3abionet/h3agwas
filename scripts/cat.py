#!/usr/bin/env python3

import sys
f = open(sys.argv[1])
col = sys.argv[2]
val = sys.argv[3]


print (f.readline)
for line in f:
    cols=line.split()
    try:
      v = float(cols[int(col)])
      cat=1 if v <=float(val) else 2
    except ValueError:
      cat=0
    cols.append(str(cat))
    print("\t".join(cols))

             
            

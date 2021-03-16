#!/usr/bin/env python3
import sys

fileI=sys.argv[1]
fileout=sys.argv[2]

readinp=open(fileI)
out=open(fileout, 'w')
outerr=open(fileout+'err', 'w')


head=readinp.readline()
out.write(head)
nbcoli=len(head.split('\t'))

for line in readinp :
  nbcol=len(line.split('\t'))
  if nbcol==nbcoli :
     out.write(line)
  else :
    outerr.write(line)

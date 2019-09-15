#!/usr/bin/env python3

# extract from bim file duplicate bim

import sys

listfilebim=sys.argv[1].split(',')
listrs=set([])
listrsdup=set([])
for File in listfilebim :
   read=open(File)
   for x in read :
      tx=x.split()
      if tx[1] in listrs :
        listrsdup.add(tx[1])
      listrs.add(tx[1])
   read.close() 

writeres=open(sys.argv[2],'w')   
writeres.writelines("\n".join(list(listrsdup)))
writeres.close()

#!/usr/bin/env python3

# extract from bim file duplicate bim

import sys

listfilebim=sys.argv[1].split(',')
listrs=set([])
listrsdup=set([])
listpos=set([])
listposdup=set([])
for File in listfilebim :
   read=open(File)
   for x in read :
      tx=x.split()
#      pos=tx[0]+" "+tx[3]
#      if pos in listpos:
#        listposdup.add(pos)
#      listpos.add(pos)
      if tx[1] in listrs :
        listrsdup.add(tx[1])
      listrs.add(tx[1])
   read.close() 

writeres=open(sys.argv[2],'w')   
writeres.writelines("\n".join(list(listrsdup)))
writeres.close()

#writeres=open(sys.argv[2]+'.dup','w')   
#writeres.writelines("\n".join(list(listposdup)))
#writeres.close()



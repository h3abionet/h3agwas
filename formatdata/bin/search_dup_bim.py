#!/usr/bin/env python3

# extract from bim file duplicate bim

import sys

listfilebim=sys.argv[1].split(',')
listrs=set([])
listinfopos=set([])
listrsdup=set([])
listpos=set([])
listposdup=set([])
writedel=open(sys.argv[3], 'w')
for File in listfilebim :
   print(File)
   read=open(File)
   for x in read :
      tx=x.split()
#      pos=tx[0]+" "+tx[3]
#      if pos in listpos:
#        listposdup.add(pos)
#      listpos.add(pos)
      if (tx[4]==tx[5]) or (len(tx[4])>1 and len(tx[5])>1):
         writedel.write(tx[0]+'\t'+tx[3]+'\t'+tx[3]+'\t'+tx[1]+'\n')
      if tx[1] in listrs :
        listrsdup.add(tx[1])
      listrs.add(tx[1])
   read.close() 

writedel.close()
writeres=open(sys.argv[2],'w')   
writeres.writelines("\n".join(list(listrsdup)))
writeres.close()

#writeres=open(sys.argv[2]+'.dup','w')   
#writeres.writelines("\n".join(list(listposdup)))
#writeres.close()



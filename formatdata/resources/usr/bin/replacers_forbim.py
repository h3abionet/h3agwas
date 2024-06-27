#!/usr/bin/env python3

# extract from bim file duplicate bim

import sys

readBim=open(sys.argv[1])
readRs=open(sys.argv[2])
writenewbim=open(sys.argv[3],'w')
ListRs=set([x.replace('\n','') for x in readRs])
listchange=set([])
for line in readBim :
    tmp=line.split()
    if len(tmp[1])>20 :
       tmp[1]=tmp[1][0:20]
    if (tmp[1] in ListRs) or tmp[1]=='.':
        if len(tmp[4])>5 :
           tmp[4]=tmp[4][0:4]
        if len(tmp[5])>5 :
           tmp[5]=tmp[5][0:4]
        newrs=tmp[0]+':'+tmp[3]+"_"+tmp[4][0:5]+"_"+tmp[5]
        Cmt=1
        while newrs in listchange :
             newrs=tmp[0]+':'+tmp[3]+"_"+tmp[4]+"_"+tmp[5]+'_'+str(Cmt)
             Cmt+=1
        tmp[1]=newrs
    listchange.add(tmp[1])
    writenewbim.write("\t".join(tmp)+'\n')

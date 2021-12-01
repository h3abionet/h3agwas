#!/usr/bin/env python3
import sys
import re
fileI=sys.argv[1]
InfoI=sys.argv[2]
PivotF=sys.argv[3]
fileout=sys.argv[4]

readpivot=open(PivotF)
headpiv=readpivot.readline().split()

DicPivot={}
for line in readpivot :
  spl=line.split()
  DicPivot[spl[0]]=line.replace('\n','')

readfile=open(InfoI)
InfoI=[x.replace('\n','') for x in readfile.readlines()]
readfile.close()
readI=open(fileI)
head=readI.readline().replace('\t\t','\t').replace('\n','').replace(" ","_").split('\t')
del head[0]
print(head)
head=head[:-2]
ncol=len(head)
head+=["P_"+x for x in InfoI]
head+=["M_"+x for x in InfoI]

writeall=open(fileout, 'w')
writeall.write("\t".join(headpiv)+"\t"+"\t".join(head)+"\n")
for l in readI :
  spl=l.replace('\n','').split()   
  rsid=spl[0]
  del spl[0]
  spl=[x for x in spl if len(x)>0]
  lines=DicPivot[rsid]+"\t"+"\t".join(spl)+'\n'
  writeall.write(lines)  
#[writeall.write(re.sub("\t$","",l)) for l in readI]
readI.close()




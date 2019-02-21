#!/usr/bin/env python3
import sys
fileI=sys.argv[1]
InfoI=sys.argv[2]
fileout=sys.argv[3]

readfile=open(InfoI)
InfoI=[x.replace('\n','') for x in readfile.readlines()]
readfile.close()
readI=open(fileI)
head=readI.readline().replace('\n','').replace(" ","_").split('\t')
head=head[:-2]
ncol=len(head)
head+=["P_"+x for x in InfoI]
head+=["M_"+x for x in InfoI]

writeall=open(fileout, 'w')
writeall.write("\t".join(head)+"\n")
[writeall.write(l) for l in readI]
readI.close()




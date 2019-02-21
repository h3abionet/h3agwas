#!/usr/bin/env python3
import sys

FileOut=sys.argv[1]
ListFile=sys.argv[2:]
NewLFile=[]
for File in ListFile :
   lire=open(File)
   tmp=lire.readline().replace('\n','').split()
   if 'CHRO' in tmp and 'POS' in tmp and 'FREQA1' in tmp and 'N' in tmp:
      NewLFile.append(File)
   lire.close()

ecrire=open(FileOut, 'w')
ecrire.write("\n".join(NewLFile))
ecrire.close()

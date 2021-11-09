#!/usr/bin/env python3
import sys

if sys.argv[1]=="NOFILE" : 
   Read=open(sys.argv[2])
   nbline=len(Read.readlines())-1
   print(nbline)
   if nbline<0 :
      print("nb line error")     
   Read.close()
   Write=open(sys.argv[3], 'w')
   #tmp=["1 0"]*nbline
   #for x in range(0,len(tmp),2):
   #   tmp[x]="1 0"
   tmp=["1"]*nbline
   Write.write("No.bed\n"+"\n".join(tmp)+"\n")
   Write.close()

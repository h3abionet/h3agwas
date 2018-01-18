#!/usr/bin/env python3

imissf  = open("$missing")

EOL=chr(10)
imiss={}
removed={}
imissf.readline()
for line in imissf:
    fields = line.strip().split()
    imiss[fields[0]+" "+fields[1]] = float(fields[5]);

genomef=open("$ibd_genome")
outf   =open("$outfname","w")
genomef.readline()
for line in genomef:
    fields = line.strip().split()
    id1 = fields[0]+" "+fields[1]
    id2 = fields[2]+" "+fields[3]
    if(float(fields[9]) > ${pi_hat}):
	if imiss[id1] >= imiss[id2] and id1 not in removed:
           outf.write("%s%s"%(id1,EOL))
           removed[id1]=1;
	elif imiss[id1]<imiss[id2] and id2 not in removed:
           outf.write("%s%s"%(id2,EOL));
           removed[id2]=1;
outf.close()

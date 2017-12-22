#!/usr/bin/env python

# Adapted from previous Perl version which had minor bugs due to 
# PLINK's visual column style -- some columns start with blank
# My Perl being rusty I rewrote to Python
 
cut_het_high=float("${params.cut_het_high}");
cut_het_low=float("${params.cut_het_low}");
cut_miss=float("${params.cut_mind}");

missf=open("$imiss");
hetfile=open("$het");
all_het=hetfile.readlines()

out=open("$outfname","w");
tab=unichr(9)
lf =unichr(10)
line=0;
for data in missf:
  data = data.strip()
  if(line>=1):
     parts_miss=data.split();
     missing=float(parts_miss[5]);
     this_het = all_het[line].strip();
     parts_het=this_het.split()
     meanHet=float(int(parts_het[4])-int(parts_het[2]))/int(parts_het[4]);
     if missing>cut_miss or meanHet>cut_het_high or meanHet<cut_het_low:
	 out.write("%s%s%s%s%s%s%s%s"%(parts_miss[0],tab,parts_miss[1],tab,missing,tab,meanHet,lf))
  line=line+1
out.close()

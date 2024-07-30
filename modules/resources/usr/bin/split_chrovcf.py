#!/usr/bin/env python3

#import bgzip
import gzip
import sys
import re



def read_header(f):
   first=""
   last =""
   contig = {}
   for line in f:
      print(line)
      #line=line.decode("utf-8")
      if "contig=" in line: break
      first = first+line
   while True:
      if "contig=" in line:
        m=re.search("contig=<ID=(.*)>",line)
        cname = m.group(1)
        contig[cname]=line
      else :
        last=last+line
      if "#CHROM" in line: break
      if line[0]!= '#' :
          print("error ")
          sys.exit(2)
      line=f.readline()#.decode("utf-8")
   return (first, contig, last)




def split(f,first,contig,last):
   chrom_files={}
   for line in f:
      data=line.split()
      the_chrom = data[0]
      if the_chrom not in chrom_files:
         raw=base+"_"+the_chrom+".vcf.gz"#"%s_%s"%(base,the_chrom),"wb")
         chrom_files[the_chrom] = gzip.open(raw ,'wb')
         chrom_files[the_chrom].write(first.encode())
         chrom_files[the_chrom].write(contig[the_chrom].encode())
         chrom_files[the_chrom].write(last.encode())
      chrom_files[the_chrom].write(line.encode())
   for f in chrom_files:
      f.close()


base  = sys.argv[1]
#f=gzip.open(input)
#input = sys.stdin
f = sys.stdin
headers = read_header(f)
split(f,*headers)

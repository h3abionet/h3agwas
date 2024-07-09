## write scott 070622
import bgzip
import gzip
import sys
import re


input = sys.stdin
base  = sys.argv[1]

def read_header(f):
   first=""
   last =""
   contig = {}
   for line in f:
      line=line.decode("utf-8")
      if "contig=" in line: break
      first = first+line
   while "contig=" in line:
      m=re.search("contig=<ID=(.*),len",line)
      cname = m.group(1)
      contig[cname]=line
      line=line.readline().decode("utf-8")
   while True:
      last=last+line
      if "#CHROM" in line: break
      line=f.readline().decode("utf-8")
   return (first, contig, last)




def split(f,first,contig,last):
   chrom_files={}
   for line in f:
      data=line.split()
      the_chrom = data[0]
      if the_chrom not in chrom_files:
         raw=open("%s_%s"%(base,the_chrom),"wb")
         chrom_files[the_chrom] = bgzip.BGZipWriter(raw)
         chrom_files[the_chrom].write(first)
         chrom_files[the_chrom].write(contig[the_chrom])
         chrom_files[the_chrom].write(last)
      chrom_files[the_chrom].write(line)
   for f in chrom_files:
      f.close()

f=gzip.open(input)
headers = read_header(f)
split(f,*headers)

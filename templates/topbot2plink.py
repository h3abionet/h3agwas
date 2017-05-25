#!/usr/bin/env python
# Takes as input
#  - A file describing an Illumica chip 
#     It should have a header line columns within the first 15 lines "Name Chr MapInfo deCODE(cM):
#     the cm is optional
#  - A file with the calls for the chip
#     There should be some header lines 
#  - the base name of the PLINK output files
#  The output l


from __future__ import print_function

import sys
import argparse
import re
def parseArguments():
    parser=argparse.ArgumentParser()
    parser.add_argument('array', type=str, metavar='samplesheet'),
    parser.add_argument('report', type=str, metavar='report',help="TOP/BOT report"),
    parser.add_argument('output', type=str, metavar='fname',help="output base"),
    args = parser.parse_args()
    return args

# auxiliary defs
chr2chr = map(str,range(0,27))
chr2chr[23]="X"
chr2chr[24]="Y"
chr2chr[25]="XY"
chr2chr[26]="MT"



def conv(x):
   try:
      num = int(x)
   except ValueError:
      if x == "X": num=23
      elif x == "Y": num=24
      elif x == "XY": num =25
      elif x == "MT": num=26
      else: num = 0
   return num

def parseArray(fname):
    f = open(fname)
    for i in range(15):
        line = f.readline()
        if "Name" in line: break
    else:
        sys.exit("Cannot find header line in "+fname)
    fields=re.split("[,\t]",line.rstrip())
    name_i = fields.index("Name")
    indices = [fields.index("Chr"),fields.index("MapInfo")]
    print(fields)
    if "deCODE(cM)" in fields:
        indices.append(fields.index("deCODE(cM)"))
    array = {}
    for line in f:
        fields=re.split("[,\t]",line.rstrip())
        curr  =[conv(fields[indices[0]]), int(fields[indices[1]])]
        if len(indices)==3:
            cm = fields[indices[2]]
            cm = 0.0 if  "NA" in cm else float(cm)
            curr.append(cm)
        array[fields[name_i]]=curr
    return array

def parseChipReport(array,fname,output):
    f = open(fname)
    for i in range(15):
        line = f.readline()
        if "SNP Name" in line: break
    else:
        sys.exit("Cannot find header line in "+fname)
    #SNP NameSample IDAllele1 - TopAllele2 - Top
    fields=re.split("[,\t]",line.rstrip())
    name_i = fields.index("SNP Name")
    samp_i = fields.index("Sample ID")
    alle_1 = fields.index("Allele1 - Top")
    alle_2 = fields.index("Allele2 - Top")
    lgenf = open ("{}.lgen".format(output),"w")
    for line in f:
        fields   = re.split("[,\t]",line.rstrip())
        snp_name = fields[name_i]
        if snp_name  not in array:
            print("Unknown SNP name in line "+line)
            continue
        a1       = fields[alle_1]
        a2       = fields[alle_2]
        lgenf.write("{}\t{}\t{}\t{}\t{}\n".format(fields[samp_i],fields[samp_i],snp_name,a1,a2))
    lgenf.close()


def outputMap(array,outname):
    entries = [[] for chrom in range(27) ]
    mapf= open("{}.map".format(outname),"w")
    i=0
    for snp in array:
        curr = array[snp]
        data = curr[1:]
        data.append(snp)
        entries[curr[0]].append(data)
        i=i+1
    for chrom in range(27):
        entries[chrom].sort()
        print(chrom,len(entries[chrom]))
        print(type(entries[chrom]))
        for [pos,cm,snp] in entries[chrom]:
            mapf.write("{}\t{}\t{}\t{}\n".format(chrom,snp,cm,pos))
    mapf.close()
            
    

args = parseArguments()

array = parseArray(args.array)
print("Done")
outputMap(array,args.output)
parseChipReport(array,args.report,args.output)



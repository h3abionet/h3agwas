#!/usr/bin/env python3
# Takes as input
#  - A file describing an Illumica chip 
#     It should have a header line columns within the first 15 lines "Name Chr MapInfo deCODE(cM):
#     the cm is optional
#  - A file with the calls for the chip
#     There should be some header lines 
#  - the base name of the PLINK output files
#  The output l

# (c) University of the Witwatersand, Johannesburg on behalf of the H3ABioinformatics Network Consortium
# 2016-2018
# Licensed under the Creative Commons Attribution 4.0 International Licence. 
# See the "LICENSE" file for details

from __future__ import print_function

import sys
import argparse
import re
import gzip
import numpy as np
from shutil import copyfile
import pandas as pd


null_values = [0,"0","",False,"false","FALSE","False"]


def parseArguments():
    parser=argparse.ArgumentParser()
    parser.add_argument('array', type=str, metavar='array description'),
    parser.add_argument('report', type=str, metavar='report',\
                        help="genotypereport"),
    parser.add_argument('samplesize', type=int, metavar='samplesize',\
                        help="how many indivs in each site, 0=all")
    parser.add_argument('idpat', type=str, metavar='idpat',help="id pattern"),
    parser.add_argument('output', type=str, metavar='fname',help="output base"),
    args = parser.parse_args()
    return args


TAB=chr(9)
EOL=chr(10)

# auxiliary defs
chr2chr = list(map(str,range(0,27)))
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
        if ",Name," in line or "Name"+TAB in line: break
    else:
        sys.exit("Cannot find header line in "+fname)
    fields=re.split("[,\t]",line.rstrip())
    name_i = fields.index("Name")
    indices = [fields.index("Chr"),fields.index("MapInfo")]
    if "deCODE(cM)" in fields:
        indices.append(fields.index("deCODE(cM)"))
    array = {}
    snp_elt=[]
    i=0
    for line in f:
        if "[Controls]" in line: break
        fields=re.split("[,\t]",line.rstrip())
        if len(indices)==3:
            cm = fields[indices[2]]
            cm = 0.0 if  "NA" in cm else float(cm)
        else:
            cm = 0.0
        snp_elt.append([conv(fields[indices[0]]), int(fields[indices[1]]), cm, fields[name_i]])
    snp_elt.sort()
    for i,content in enumerate(snp_elt):
        array[content[-1]]=i
    return snp_elt, array

def generate_line(pedf,old_sample_id,output):
    pedf.write(TAB.join((old_sample_id,"0","0","0","0")))
    pedf.write(TAB)
    pedf.write(TAB.join(output)+EOL)
    pass

def getReportIndices(line):
    #SNP NameSample IDAllele1 - TopAllele2 - Top
    fields=re.split("[,\t]",line.rstrip())
    #fields = line.rstrip().split(",")
    if "Sample Name" in fields:
        sample = "Sample Name"
    elif "Sample ID" in fields:
        sample = "Sample ID"
    else:
        sys.exit("Can't find sample id as a field in <%s>"%line)
    if "Allele1 - Top" in fields:
        allele1 = "Allele1 - Top"
        allele2 = "Allele2 - Top"
    elif "Allele1 - Forward" in fields:
        allele1, allele2 = "Allele1 - Forward","Allele2 - Forward"
    else:
        sys.exit("Can't find field labels to mark alleles in <%s>"%line)
    name_i = fields.index("SNP Name")
    samp_i = fields.index(sample)
    alle_1 = fields.index(allele1)
    alle_2 = fields.index(allele2)
    return name_i, samp_i, alle_1, alle_2

def getID(idreg,sample_id):
    m = idreg.match(sample_id)
    if sample_id in replicates:
        rep = "_replicate_"+replicates[sample_id]
    else:
        rep = ""
    if m:
       if len(m.groups())==2:
           return (m.group(1)+rep, m.group(2))
       elif len(m.groups())==1:
           return (m.group(1)+rep,m.group(1))
       else:
           sys.exit("Pattern <%s> has wrong number of groups"%args.idpat)
    else:
       sys.exit("Sample ID <%s> cannot be parsed by <%s>"%(sample_id,args.idpat))
    return sample_id

def checkgzipfile(fname) :
    GZIP_MAGIC_NUMBER = "1f8b"
    f = open(fname)
    isGz=f.read(2).encode("utf-8").hex() == GZIP_MAGIC_NUMBER
    f.close()
    return isGz


def parseChipReport(snp_elt,array,fname,output):
    # how many lines do we need to skip
    # Looks like 10, but let's be sure
    idreg = re.compile(args.idpat)
    if checkgzipfile(fname) :
      f = gzip.open(fname,"rt")
    else :
      f = open(fname)
    head=0
    for line in f:
        head=head+1
        if "[Data]" in line: break
    name_i, samp_i, alle_1, alle_2 = getReportIndices(f.readline())
    pedf = open ("{}.ped".format(output),"w")
    old_sample_id="xxx" # ("xxx","")
    num=0
    output = np.empty([len(snp_elt)*2],dtype='U1')
    output.fill("0")
    for line in f:
        #fields = line.rstrip().split(",")
        fields = re.split("[,\t]",line.rstrip())
        sample_id = fields[samp_i].replace(' ','') 
        snp_name = fields[name_i]
        if snp_name  not in array:
            sys.exit("Unknown SNP name in line "+line)
        a1       = fields[alle_1]
        a2       = fields[alle_2]
        if a1 == "-": a1="0"
        if a2 == "-": a2="0"
        if sample_id != old_sample_id:
            if num > args.samplesize >  0:
                pedf.close()
                return
            if old_sample_id!="xxx":
                generate_line(pedf,old_sample_id,output)
            output.fill("0")
            old_sample_id = sample_id
            num=num+1
        ind = array.get(snp_name)
        output[2*ind]=a1
        output[2*ind+1]=a2
    generate_line(pedf,old_sample_id,output)
    pedf.close()


def outputMap(snp_elt,array,outname):
    mapf= open("{}.map".format(outname),"w")
    for [chrom,pos,cm,snp] in snp_elt:
        mapf.write("{}{}{}{}{}{}{}{}".format(chrom,TAB,snp,TAB,cm,TAB,pos,EOL))
    mapf.close()
            
if len(sys.argv) == 1:
   sys.argv=["topbot2plink.py","$array","$report","$samplesize", "$idpat", "$mask", "$replicates", "$output"]
    


args = parseArguments()
if args.idpat in null_values:
    args.idpat=("(.*)")



snp_elt, array = parseArray(args.array)
parseChipReport(snp_elt,array,args.report,args.output)

outputMap(snp_elt,array,args.output)


#copyFam(args.fam,args.output)


#!/usr/bin/env python

# Converts from CRLMM format into PLINK format
# Takes two arguments -- name of CRLMM file and the base of the PLINK file
#    
# Tested under both Python 2.7 and 3.5.2
# Python 2 seems significantly faster
# 
# Scott Hazelhurst on behalf of the H3ABioinformatics Network Consortium
# December 2016
# (c) Released under GPL v.2


from __future__ import print_function

import sys
import os
import re

# check if we are being called from the command line or as a template
# from Nextflow. If  as a template from Nextflow, then the nextflow process
# must define CRLMMCALL and OUTBASE variables
# NB: only for a Nextflow template -- you could also directly call from a Nextflow
# script if this script is on the path -- in which case the parameters should be passed
if len(sys.argv)==0:
    sys.argv = [sys.argv[0],"$CRLMMCALL","$MANIFEST","$OUTBASE"]

max_buffer_size = 1<<20


# This is the PLINK magic number and also the SNP-major mode byte
magic_number_mode=[0b01101100,0b00011011,0b00000001]

# Conversion from code that CRLMM uses for homo/heterozygosity to what PLINK used
crlmm2bed = { "0":2, "1":0, "2":1, "3":3 }


def produceFam(fin,base):
    line = fin.readline().replace('"','').rstrip().split()
    outf = open("%s.fam"%base,"w")
    for pid in line:
        outf.write("%s\t%s\t0\t0\t0\t0\n"%(pid,pid))
    outf.close()

def alleles(data):
    m = re.search(".(.*)/(.*).",data)
    A = m.group(1)
    B = m.group(2)
    return (A,B)


def getManifest(mfile):
    snpdet={}
    with open(mfile) as f:
        f.readline()
        for line in f:
            data = line.split()
            snpdet[data[1]]=(alleles(data[3]),data[9],data[10])
    return snpdet

def checkPlainSNV(data):
       ordinary_snp=True
       for d in data[1:]:
           if d not in ["0","1","2","3"]: ordinary_snp = False
       return ordinary_snp

# Given a file in CRLMM format, returns the corresponding PLINK-encoded
# lines in SNP-major format
def bufferise(fin,bimf):
    n=0
    for snp in fin:
        n=n+1
        bval=0 # byte we are building up
        i=0      # which 2-bits we're on
        data = snp.rstrip().split()
        # we can only handle biallelic SNPs
        if not checkPlainSNV(data): continue
        snp = data[0].replace('"','')
        try:
            ((A,B),chrom,pos)=snpdet[snp]
            bimf.write("%s\t%s\t0\t%s\t%s\t%s\n"%(chrom,snp,pos,A,B))
        except KeyError:
            print("Warning no manifest for SNP %s "%snp)
            bimf.write("%s\t%s\t0\t%s\t%s\t%s\n"%(0,snp,0,0,0))
        # go through each person's allele for the current SNP
        for person in data[1:]:
            bval = (bval << 2)|crlmm2bed[person]
            i=i+1
            if  i== 4:
                yield bval
                bval=i=0
        if i!=0: yield bval


def produceBedBim(fin,base):
    bedf = open("%s.bed"%base,"wb")
    bimf = open("%s.bim"%base,"w")
    bedf.write(bytearray(magic_number_mode))
    buffer=[]
    # get each byte in turn from input  and store in 'buffer'
    for bval in bufferise(fin,bimf):
         buffer.append(bval)
         # at regular intervals sump buffer to disk
         if len(buffer)==max_buffer_size:
              bedf.write(bytearray(buffer))
              buffer=[]
    # anything at end needs to be dumped to disk
    if len(buffer)>0:
         bedf.write(bytearray(buffer))
    bedf.close()
    bimf.close()

if __name__ == "__main__":                 
    fin = open(sys.argv[1])
    snpdet = getManifest(sys.argv[2])
    base = sys.argv[3]
    produceFam(fin,base)
    produceBedBim(fin,base)









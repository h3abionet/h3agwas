#!/usr/bin/env python3
author__ = 'glebkichaev'

import sys
import numpy as np
from bisect import *
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-c", "--chr", dest="chr_name")
parser.add_option("-p", "--pos", dest="pos_name")
parser.add_option("-i", "--input", dest="input_file")
parser.add_option("-o", "--out", dest="output_file")
parser.add_option("-l", "--locus", dest="locus_name")

(options, args) = parser.parse_args()

chr_name=options.chr_name
pos_name=options.pos_name
input_file=options.input_file
output_file=options.output_file
locus_name=options.locus_name

usage = \
"""Usage:
    -i (--input) file that contains list of paths for annotation .bed files
    -l (--locus) file with chromosome and positions of snps in region to be annotated 
    -o (--out) name of filename to output
    -c (--chr) identifier of chromosome field in the header
    -p (--pos) identifier of position fieild in the header (hg 19)""" 


if chr_name == None or pos_name == None or input_file == None or output_file == None or locus_name==None:
    sys.exit(usage)


annotation_file = open(input_file)
locus_file = open(locus_name)

header = locus_file.readline()
header = header.strip().split()
chr_ind = header.index(chr_name)
pos_ind = header.index(pos_name)
chr = []
positions = []
for snps in locus_file:
    temp = snps.strip().split()
    positions.append(int(float(temp[pos_ind])))
    chr.append(temp[chr_ind])

all_annotations = []
all_names = []
for annotations in annotation_file:
    print(annotations)
    current_annotation = open(annotations.strip())
    name = annotations.strip().split("/")[-1]
    chr_locs = []
    left_int = []
    right_int = []
    for lines in current_annotation:
        site = lines.strip().split()
        if(site[0] == chr[1]):
            chr_locs.append(site)
            left_int.append(int(site[1]))
            right_int.append(int(site[2]))
            
    locus_annotation =[]
    
    tot_coordinates= len(left_int)
    for snps in positions:
        left_st = max(bisect_left(left_int, snps)-1,0)
        right_end = min(bisect_right(right_int, snps)+1,tot_coordinates)
        check_interval = [int(locs[1]) <= int(snps) <= int(locs[2]) for locs in chr_locs[left_st:right_end]] 
        if any(check_interval) == True:
            annotate_snp = 1
        else:
            annotate_snp = 0
        locus_annotation.append(annotate_snp)
    if sum(locus_annotation)>0 :
       all_annotations.append(locus_annotation)
       all_names.append(name)

out_file = open(output_file, "w")
outline = " ".join(all_names)
out_file.write(outline+"\n")
all_annotations_transpose = list(map(list, zip(*all_annotations)))
for lines in all_annotations_transpose:
    outline = " ".join(map(str,lines))
    out_file.write(outline+"\n")
out_file.close()


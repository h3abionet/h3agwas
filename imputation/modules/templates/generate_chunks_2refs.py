#!/usr/bin/env python2.7

import argparse,sys
import time

parser = argparse.ArgumentParser()
parser.add_argument("--legend_file_list", help="Searated by comma")
parser.add_argument("--outfile", help="")
parser.add_argument("--chunk_size", help="")
parser.add_argument("--chrm", help="chromosome")
args = parser.parse_args()


def chunk_split(legend_file_list, outfile, chunk_size='', chrm=''):
    '''
    Return: chunk files in the output folder
    '''
    pop_samples={}
    print 'Generating chunk files from reference data ...'
    print "Reading file(s) "+legend_file_list
    print "Writing file "+outfile
    POS = set()
    for legend_file in legend_file_list.split(','):
        for line in open(legend_file):
            if "id" not in line and 'position' not in line:
                line = line.strip().split()
                try:
                    pos = int(line[1])
                    POS.add(pos)
                except:
                    pass
    max_POS = max(POS)
    chunk_size = int(chunk_size)
    out=open(outfile,"wt")
    # out.writelines("start" +" "+"end"+"\n")
    for pos in list(range(1, max_POS+chunk_size, chunk_size)):
        start_ = pos
        end_ = start_ + chunk_size - 1
        out.writelines(str(chrm)+","+str(start_) +","+str(end_)+"\n")
    out.close()

chunk_split(args.legend_file_list, args.outfile, args.chunk_size, args.chrm)
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

def parseArguments():
    parser=argparse.ArgumentParser()
    parser.add_argument('array', type=str, metavar='samplesheet'),
    parser.add_argument('report', type=str, metavar='report',help="TOP/BOT report"),
    parser.add_argument('output', type=str, metavar='fname',help="output base"),
    args = parser.parse_args()
    return args


def parseArray(fname):
    f = open(fname)
    for i in range(15):
        line = f.readline()
        if "Name," in line: break
    else:
        sys.exit("Cannot find header line in "+fname)
    fields=line.split(",")
    indices  = map(lambda field_name : fields.index(field_name), ["Name","Chr","MapInfo","deCODE(cm)"])
    if indices[-1]=="-1": indices.pop
    indices = slice(indices)
    array={}
    for line in f:
        fields=line.split(",")[indices]
        array[fields[0]]=fields[1:]
    return array


parseArray(array)

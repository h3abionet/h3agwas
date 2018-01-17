#!/usr/bin/env python3

import argparse
import openpyxl
import sys
import re
import shutil
import os
# Expects four arguments
# -- xlsx file with the data
# -- the original fam file
# -- the output name (just the base -- e.g. chip, not chip.fam)
# -- abbrev?  If "0" then the IDs are used as provided in the samplesheet
#                 If "1", we remove everything up to the last underscore in the name


def parseArguments():
    parser=argparse.ArgumentParser()
    parser.add_argument('samplesheet', type=str, metavar='samplesheet'),
    parser.add_argument('fam', type=str, metavar='report',help="fam file"),
    parser.add_argument('output', type=str, metavar='fname',help="output base"),
    parser.add_argument('abbrev', type=str, metavar='abbrev',help="abbreviate IDs"),
    args = parser.parse_args()
    return args


# If called with no parameters, we assume using Nexflow's template mechanism and the parameters
# are substituted in with fixed names
if len(sys.argv)<=1:
    sys.argv = ["sheet2fam.py","$samplesheet","$fam", "$output","$abbrev"]


# we avoid backslashes
TAB=chr(9)
EOL=chr(10)

def getHeading(allrows):
    heading = list(map(lambda x: x.value, allrows.__next__()))
    col_id=heading.index("Institute Sample Label")
    col_sex=heading.index("Manifest Gender")
    col_batch=heading.index("Batch Comment")
    return (col_id, col_sex, col_batch)

def sex_code(x):
    if x == "Male":
        return 1
    elif x == "Female":
        return 2
    else:
        return 0


def parseSheet(allrows):
    [col_id, col_sex, col_batch]  = getHeading(allrows)
    indivs = {}
    for row in allrows:
       sample_id   = row[col_id].value
       if args.abbrev:
            m = re.search(".*_(.*)",sample_id)
            if m:
                sample_id=m.group(1)
            else:
                print("Sample ID <%s> cannot be abbrev"%sample_id)
       sample_sex  = row[col_sex].value
       batch       = row[col_batch].value
       m=re.search(".* (.+)",batch)
       batch = m.group(1)
       indivs[sample_id] =  [sample_sex,batch]
    return indivs

def produceFam(indivs,origfam):
    g = open("{}.fam".format(args.output),"w")
    for sample_id in origfam:
        [sample_sex,batch] = indivs[sample_id]
        g.write("{}{}{}{}0{}0{}{}{}{}{}".format(sample_id,TAB,sample_id,TAB,TAB,TAB,sex_code(sample_sex),TAB,batch,EOL))
    g.close()

def getFam(fam):
    """ get the list of IDs from the orignal fam file -- we need to know the correct order of the fam file """
    f = open(fam)
    origfam = []
    for line in f:
        data = line.split()
        fid=data[0]
        iid=data[1]
        origfam.append(fid)
    return origfam

args        = parseArguments()
shutil.copyfile(args.fam,"oldfam.fam")
xlsxf       = openpyxl.load_workbook(args.samplesheet)
sheet_names = xlsxf.get_sheet_names()
allrows     = xlsxf[sheet_names[0]].rows
indivs      =  parseSheet(allrows)
origfam     = getFam(args.fam)
produceFam(indivs,origfam)
xlsxf.close()



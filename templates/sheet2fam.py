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
    parser.add_argument('batch_col', type=str, metavar='batch',help="batch column"),
    parser.add_argument('output', type=str, metavar='fname',help="output base"),
    parser.add_argument('idpat', type=str, metavar='idpat',help="abbreviate IDs"),
    args = parser.parse_args()
    return args


# If called with no parameters, we assume using Nexflow's template mechanism and the parameters
# are substituted in with fixed names
if len(sys.argv)<=1:
    sys.argv = ["sheet2fam.py","$samplesheet","$fam", "$batch_col", "$output","$idpat"]


# we avoid backslashes
TAB=chr(9)
EOL=chr(10)

def getHeading(allrows):
    heading = list(map(lambda x: x.value, allrows.__next__()))
    col_id=heading.index("Institute Sample Label")
    col_sex=heading.index("Manifest Gender")
    if args.batch_col in ["0",0,False,"false",None,""]:
        col_batch=-1
    else:
        col_batch=heading.index(args.batch_col)
    return (col_id, col_sex, col_batch)

def sex_code(x):
    if x == "Male":
        return "1"
    elif x == "Female":
        return "2"
    else:
        return "0"


def parseSheet(allrows):
    [col_id, col_sex, col_batch]  = getHeading(allrows)
    indivs = {}
    batch="-9"
    for row in allrows:
       sample_id   = row[col_id].value
       if args.idpat not in [0, "0", False, ""]:
            m = re.search(args.idpat,sample_id)
            if m:
                sample_id=m.groups()
                if len(sample_id)==1: sample_id=sample_id+sample_id
            else:
                print("Sample ID <%s> cannot be abbrev"%sample_id)
       else:
           sample_id=(sample_id,sample_id)
       sample_sex  = row[col_sex].value
       if col_batch>0:
           batch       = row[col_batch].value
           m=re.search("Batch (.+)",batch)
           if m:
               batch = m.group(1)
       indivs[sample_id] =  [sample_sex,batch]
    return indivs

def produceFam(indivs,origfam):
    g = open("{}.fam".format(args.output),"w")
    for sample_id in origfam:
        [sample_sex,batch] = indivs[sample_id]
        data = TAB.join(sample_id+("0","0",sex_code(sample_sex),batch))+EOL
        g.write(data)
    g.close()

def getFam(fam):
    """ get the list of IDs from the orignal fam file -- we need to know the correct order of the fam file """
    f = open(fam)
    origfam = []
    for line in f:
        data = line.split()
        fid=data[0]
        iid=data[1]
        origfam.append((fid,iid))
    return origfam

args        = parseArguments()
if args.samplesheet in [0,"","0",False]:
    shutil.copyfile(args.fam,"${output}.fam")
    sys.exit(0)

xlsxf       = openpyxl.load_workbook(args.samplesheet)
sheet_names = xlsxf.get_sheet_names()
allrows     = xlsxf[sheet_names[0]].rows
indivs      =  parseSheet(allrows)
shutil.copyfile(args.fam,"oldfam.fam")

origfam     = getFam(args.fam)
produceFam(indivs,origfam)
xlsxf.close()



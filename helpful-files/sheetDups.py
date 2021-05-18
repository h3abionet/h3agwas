#!/usr/bin/env python3

import argparse
import openpyxl
import sys
import re
import shutil
import os
import pandas as pd
# Expects four arguments
# -- xlsx file with the data
# -- the original fam file
# -- the output name (just the base -- e.g. chip, not chip.fam)
# -- abbrev?  If "0" then the IDs are used as provided in the samplesheet
#                 If "1", we remove everything up to the last underscore in the name



def parseArguments():
    parser=argparse.ArgumentParser()
    parser.add_argument('--batch-col', dest="batch_col", type=str, metavar='batch',\
                        help="batch column"),
    parser.add_argument('samplesheet', type=str, metavar='samplesheet',nargs='+'),
    parser.add_argument('--bad-batch', dest="bad_batch", default=False)
    parser.add_argument('--phe', dest="phe", type=str, metavar='phe',\
                        help="phe"),
    parser.add_argument('--output', dest="output", type=str, metavar='fname',\
                        help="output base"),
    parser.add_argument('--idpat', dest="idpat", type=str, metavar='idpat',\
                        help="abbreviate IDs"),
    parser.add_argument('--out_idpat', dest="out_idpat", type=str, metavar='out_idpat',\
                        help="abbreviate IDs"),
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
    call    = heading.index("Call_Rate")
    if "Institute Plate Label" in heading:
        label = "Institute Plate Label"
    else:
        label = "Sample Plate"
    plate   = heading.index(label)
    well    = heading.index("Well")
    # Now we get the sex -- used to be labelled "Gender" but we asked Illumina to change but old sheets will
    # have this
    if "Manifest Gender" in heading:
        sex_colname = "Manifest Gender"
        sex_geno    = "Genotype Gender"
    elif "Manifest Sex" in heading:
        sex_colname = "Manifest Sex"
        sex_geno    = "Genotyped Sex"
    else:
        sys.exit("Can't find manifest sex column in sample sheet <%s>"%args.samplesheet)
    col_sex=heading.index(sex_colname)
    sex_geno_col = heading.index(sex_geno)
    if args.batch_col in ["0",0,False,"false",None,""]:
        col_batch=-1
    else:
        col_batch=heading.index(args.batch_col)
    return (col_id, col_sex, sex_geno_col, col_batch,call,plate,well)

def sex_code(x):
    if x == "Male":
        return "1"
    elif x == "Female":
        return "2"
    else:
        return "0"


def parseSheet(allrows,indivs,sofar,problems):
    [col_id, col_sex, col_geno_sex, col_batch,call,plate,well]  = getHeading(allrows)
    batch="-9"
    for row in allrows:
       sample_id   = raw_id = row[col_id].value.replace(" ","")
       if sample_id in exclusions : pass
       if col_batch>0:
           batch       = row[col_batch].value
           m=re.search("Batch (.+)",batch)
           if m:
               batch = int(m.group(1))
       m = re.search(args.out_idpat,sample_id)
       if m:
           outid  = m.group(1)
       else:
           sys.exit("Out ID PAT does not match on "+sample_id)
       if args.idpat not in [0, "0", False, ""]:
            m = re.search(args.idpat,sample_id)
            if m:
                sample_id=m.groups()
                if len(sample_id)==1: sample_id=sample_id+sample_id
            else:
                print("Sample ID <%s> cannot be abbrev"%sample_id)
       else:
           sample_id=(sample_id,sample_id)
       data = [sample_id,raw_id, row[col_sex].value, row[col_geno_sex].value,\
               batch, row[call].value, row[plate].value, row[well].value]
       if outid in sofar:
           if outid in problems:
               problems[outid] = problems[outid]+[data]
           else:
               problems[outid] = [sofar[outid],data]
       sofar[outid]=data
       indivs[outid]=1






def createReplicates(fname,fd,problems):
    g = open("%s.rep"%fname,"w")
    h = open("%s.err"%fname,"w")
    m = open("%s.miss"%fname,"w")
    allf = open("%s.all"%fname,"w")
    for k in sorted(problems.keys()):
        try:
            fam_sex = str(fd.loc[k]['sex'])
        except KeyError:
            fam_sex = 0
        best=-1
        rate=0
        for i, v in enumerate(problems[k]):
            sex_ok = (fam_sex==sex_code(v[3])) or  (sex_code(v[3])==sex_code(v[2])) and ("0"==fam_sex)
            if (v[4]>rate) and sex_ok :
                rate=v[5]
                best=i
        dup=edup=1
        ok=False
        errs=[]
        allf.write('%s '%k)
        allf.write(' '.join(map(str,sofar[k][1:]))+EOL)
        for i, v in enumerate(problems[k]):
            if i==best: 
                ok=True
                continue
            sex_ok = (fam_sex==sex_code(v[3])) or  (sex_code(v[3])==sex_code(v[2])) and ("0"==fam_sex)
            if sex_ok:
                g.write("%s%s%s%s%s%s%d%s"%(v[0][0],TAB,v[0][1],TAB,v[1],TAB,dup,EOL))
                allf.write("%s_replicate_%d "%(k,dup)+" ".join(map(str,v[1:]))+EOL)
                dup=dup+1
            else:
                errs.append("%s%s%s%s%s%s%d%s"%(v[0][0],TAB,v[0][1],TAB,v[1],TAB,edup,EOL))
                edup=edup+1
        if not ok:
            m.write(errs[0])
            del errs[0]
        h.writelines(errs)
    for x in exclusions:
        h.write("%s\t%s\n"%(x,x))
    g.close()
    h.close()
    m.close()
    allf.close()
        


args        = parseArguments()


indivs   = {}
problems = {}
sofar    = {}

if args.bad_batch:
    f = open(args.bad_batch)
    exclusions = set(map(lambda x:x.strip(), f.readlines()))
else:
    exclusions=[]

for name in args.samplesheet:
    xlsxf       = openpyxl.load_workbook(name)
    sheet_names = xlsxf.sheetnames
    allrows     = xlsxf[sheet_names[0]].rows
    parseSheet(allrows,indivs,sofar,problems)
xlsxf.close()

if args.phe:
    fd = pd.read_csv(args.phe,delim_whitespace=True,index_col="FID")
else:
    fd = pd.DataFrame.from_dict(indivs,orient='index',columns=["sex"])



if len(problems)>0:
    print("ID,Sample Label,Manifest Sex,Genotype Sex,Batch,Call Rate,Plate,Well")
    for k in sorted(problems.keys()):
        #print(problems[k])
        pass
    if args.output:
        createReplicates(args.output,fd,problems)




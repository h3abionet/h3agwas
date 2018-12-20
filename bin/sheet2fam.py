#!/usr/bin/env python3
# (c) University of the Witwatersand, Johannesburg on behalf of the H3ABioinformatics Network Consortium
# 2016-2018
# Licensed under the Creative Commons Attribution 4.0 International Licence. 
# See the "LICENSE" file for details

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
# optionally files with list of replicat
#    replicates -- they have the sample label and a number for the replicate


null_values = [0,"","0",False,"False","false","FALSE"]

def parseArguments():
    parser=argparse.ArgumentParser()
    parser.add_argument('--idpat', dest='idpat',default=False, metavar='idpat',help="1st pattern match"),
    parser.add_argument('--newpat', dest='newpat',default=False, metavar='newpat',help="2nd pattern match"),
    parser.add_argument('samplesheet', type=str, metavar='samplesheet'),
    parser.add_argument('fam', type=str, metavar='report',help="fam file"),
    parser.add_argument('--sheet-columns', type=str, metavar='fname',default="0",\
                        dest='sheet',help="file with mappings for columns"),
    parser.add_argument('--replicates', metavar='replicates',default="0",\
                        dest='replicates',help="file with mappings for columns"),
    parser.add_argument('--mask', metavar='mask',default="0",\
                        dest='mask',help="file with mappings for columns"),
    parser.add_argument('output', type=str, metavar='fname',help="output base"),
    args = parser.parse_args()
    return args




# we avoid backslashes
TAB=chr(9)
EOL=chr(10)

legal_columns = ['Institute Sample Label','Sample Plate','Well','Manifest Sex',"Batch Comment"]
column_index  = ['sample_label','plate','well','sex','batch']


def getSheetColumnMaps(fname):
    column = dict(zip(column_index,legal_columns))
    if fname in null_values: return column
    f = open(fname)
    for line in f:
        try:
            (var,val) = line.split('=')
        except ValueError:
            sys.exit("Problem with file <%s> line <%s>: no valid definition "%(fname,line))
        var = var.strip()
        if var not in column_index:
            sys.exit("Failed reading in sheet config in file <%s> : don't know variable <%s>"%(fname,var))
        val = val.strip()
        column[var]=val if val !=0 else -1
    return column
    



def extractCol(heading,col):
    try:
        if ".xls" in args.samplesheet:
            result = heading.index(column[col])
        else:
            result = column[col]
    except (ValueError, KeyError) as e:
        sys.exit("You have used a label <%s> but not in <%s>"%(col,str(heading)))
    return result

def getHeading(allrows):
    # refactor to use pandas instead of openpyxl
    if ".xls" in args.samplesheet:
        heading = list(map(lambda x: x.value, allrows.__next__()))
    else:
        heading = allrows
    col_id    =  extractCol(heading,'sample_label')
    col_plate =  extractCol(heading,'plate')
    col_well  =  extractCol(heading,'well')
    col_sex   =  extractCol(heading,'sex')
    col_batch =  extractCol(heading,'batch')
    return (col_id, col_plate, col_well, col_sex, col_batch)

def sex_code(x):
    if x == "Male":
        return "1"
    elif x == "Female":
        return "2"
    else:
        return "0"


def getID(pattern, raw_id):
   if pattern not in null_values:
      m = re.search(pattern,raw_id)
      if m:
         sample_id=m.groups()
         if len(sample_id)==1: sample_id=sample_id+sample_id
      else:
            sys.exit("Sample ID <%s> cannot be abbreviated by <%s>"%(raw_id,pattern))
   else:
       sample_id=(raw_id,raw_id)
   fid, iid = sample_id
   return (fid.strip(),iid.strip())



def getIndivHash(fname):
    """ Open a file and build a dict of individuals -- assume that first column is the individual ID
        and the last column is a label we might be interested in -- e.g. replicate number  """
    indivs = {}
    if fname in null_values: return indivs
    f = open(fname)
    for line in f:
        if line.strip() == "0" : return {}
        data = line.rstrip().split()
        indivs[(data[0],data[1])]=data[-1]
    return indivs


def parseSheet(allrows):
    [col_id, col_plate,col_well, col_sex, col_batch]  = getHeading(allrows)
    sofar = {}
    problems = {}
    indivs = {}
    batch="-9"
    def getVal(row,col):
       try:
          if ".xls" in args.samplesheet:
             the_val=row[col].value
             if the_val is None: # if the cell is empy
                 the_val = ""
             elif type(the_val) in [int,float]:
                 the_val=str(the_val)
             return the_val.replace(" ","")
          else:
             return row[col].replace(" ","")
       except KeyError as e:
          print(EOL+"<%s> is not a column of the sample sheet"%col+EOL)
    for row in allrows:
       raw_id = getVal(row,col_id)
       (fid,iid)= getID(args.idpat, raw_id)
       if (fid,iid) in masks: continue
       if args.newpat not in null_values and (fid != getVal(row,col_plate) or iid != getVal(row,col_well)):
           print("mismatch in this row <%s %s %s %s %s>"%(raw_id,fid,iid,getVal(row,col_plate),getVal(row,col_well)))
       (real_fid, real_iid) = getID(args.newpat, raw_id)
       print(raw_id,fid,iid,real_fid,real_iid)
       if (fid,iid)  in replicates:
           real_fid = real_fid + "_replicate_" + replicates[(fid,iid)]
       if real_fid in sofar:
           if real_fid in problems:
               problems[real_fid] = problems[real_fid]+","+raw_id
           else:
               problems[real_fid] = sofar[real_fid]+","+raw_id
       sofar[real_fid]=raw_id
       sample_sex  = getVal(row,col_sex)
       if col_batch not in null_values:
           batch       = getVal(row,col_batch)
           m=re.search("Batch (.+)",batch)
           if m:
               batch = m.group(1)
           else:
               batch=batch.replace(" ","")
       indivs[(fid,iid)] =  [real_fid,real_iid,sample_sex,batch]
    if len(problems)>0:
        print(EOL+EOL+"==============================================="+EOL+EOL)
        print("The following IDs are duplicated in the sample sheet --- a problem")
        print("If a genuine replicate, one should have the ID modifed (dd '-_replicate_')")
        for idn in problems:
            print("    ",idn,problems[idn])
        print("""" 
                The IDs above duplicated in the sample sheet --- a problem.
                Either set the "replicates" option or manually add '-_replicateN_') to the sheet ")
              """)
    return indivs, problems

def produceFam(indivs,problems,origfam):
    g = open("{}.fam".format(args.output),"w")
    for sample_id in origfam:
        try:
            (fid,iid) = getID(args.idpat, sample_id)
            ofid,oiid = fid,iid
            [fid,real_id,sample_sex,batch] = indivs[(fid,iid)]
            if (ofid,oiid) != (fid,real_id):
               print("<%s,%s> --->  <%s,%s>"%(ofid,oiid,fid,real_id))
            if real_id in problems:
                print("The ID <%s> with fid, iid, real_id <%s> <%s> <%s> is a duplicate"%(sample_id,fid,iid,real_id))
                sys.exit(124)
        except KeyError:
            m = re.search("(.)-_replicate_",sample_id)
            if m and m.group(1) in indivs:
                [fid, real_id,sample_sex,batch] = indivs[m.group(1)]
                continue
            print(EOL+EOL+EOL+"======================================================"+EOL+EOL+EOL)
            print("Your data from the genotyping report contains the ID {} -- but this is not in the sample sheet <{}>".\
                  format(sample_id,args.samplesheet))
            print(EOL+EOL+EOL+"======================================================"+EOL+EOL+EOL)
            sys.exit(125)
        data = TAB.join((fid,real_id,"0","0",sex_code(sample_sex),batch))+EOL
        g.write(data)
    g.close()

def getFam(fam):
    """ get the list of IDs from the orignal fam file -- we need to know the correct order of the fam file """
    f = open(fam)
    origfam = []
    for line in f:
        data = line.split()
        fid=data[0]
        origfam.append(fid)
    return origfam


def getCSV(fname):
    # is there a [Data]
    f=open(fname)
    skip=0
    for line in f:
        skip=skip+1
        if "[Data]" in line:
            break
    else: # Hope that the first line is the header line
        skip=0
    df = pd.read_csv(fname,dtype={'Sample Group':str,'Gender':str},skiprows=skip,delimiter=",")
    df.fillna('',inplace=True)
    return map(lambda x:x[1],df.iterrows())

args        = parseArguments()
if args.samplesheet in null_values:
    shutil.copyfile(args.fam,"${output}.fam")
    sys.exit(0)

column = getSheetColumnMaps(args.sheet) # How are the columns in the sample sheet labelled?

replicates = {} if args.replicates in null_values else getIndivHash(args.replicates)
masks      = {} if args.mask in null_values else getIndivHash(args.mask)


origfam     = getFam(args.fam)

if ".xls" in args.samplesheet:
    xlsxf       = openpyxl.load_workbook(args.samplesheet)
    sheet_names = xlsxf.sheetnames
    allrows     = xlsxf[sheet_names[0]].rows
elif ".csv" in args.samplesheet:
    allrows = getCSV(args.samplesheet)

indivs, problems      =  parseSheet(allrows)
shutil.copyfile(args.fam,"oldfam.fam")


produceFam(indivs,problems,origfam)




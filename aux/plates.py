# Scott Hazelhurst, University of the Witwatersrand, Johannesburg  (C) 2018 on behalf of H3ABioNet

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import argparse
import re

# if sample-sheet given, extract out plate, well from sample sheet. 
# else, assume that the IDs have the sample sheet embedded

def parseArguments():
    parser=argparse.ArgumentParser()
    parser.add_argument('--sample-sheet', dest="sample_sheet", \
                        type=str, metavar='batch',\
                        help="samplesheet"),
    parser.add_argument('target', help="file with IDs"),
    parser.add_argument('out',  type=str, help="out"),
    args = parser.parse_args()
    return args

args=parseArguments()

sex_f = open(args.target)
the_ids=[]
for line in sex_f:
    data = line.strip().split()
    the_ids.append(data[0])


if args.sample_sheet:
    batches = pd.read_excel(args.samplesheet)
    batches['ID'] = batches['Institute Sample Label'].apply(lambda x:x[18:])
    batches['row']=batches['Well'].apply(lambda x:ord(x[0])-ord('A'))
    batches['col']=batches['Well'].apply(lambda x:int(x[1:4]))
    problems = batches[batches['ID'].isin(the_ids)]
else:
    plist=[]
    for the_id in the_ids:
        m=re.search("(.*)_(.)(..)",the_id)
        plate=m.group(1)
        row=ord(m.group(2))-ord('A')
        col=int(m.group(3))
        plist.append((plate,row,col))
    problems=pd.DataFrame.from_records(plist,\
                                     columns=["Institute Plate Label", "row", "col"])
            

plates=problems.groupby('Institute Plate Label')

for plate_name,plate in plates:
    plt.xlim(0,13)
    plt.ylim(-1,8)
    plt.xticks(range(13),["",1,2,3,4,5,6,7,8,9,10,11,12])
    plt.yticks(range(8),["A","B","C","D","E","F","G","H"])
    plt.scatter(plate['col'],plate['row'],color='blue')
    plt.savefig("%s/%s.png"%(args.out,plate_name))
    plt.close()

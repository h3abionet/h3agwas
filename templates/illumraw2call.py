#! /usr/bin/env python

from __future__ import print_function

import argparse
import os
import sys

EOL=unichr(10)

def parseArguments():
    parser=argparse.ArgumentParser()
    parser.add_argument('sample', type=str, metavar='samplesheet'),
    parser.add_argument('idat', type=str, metavar='IDATDIR',help="directory where IDAT files can be found"),
    parser.add_argument('--sample-sample-column',dest="sample_col",type=int,default=0,\
                        help="col in sample file where sample ID can be found (number cols from 0)")
    parser.add_argument("--sentrix-barcode-column",dest="barcode_col",type=int,default=3,\
                        help="col in sample file where Sentrix barccode found (number cols from 0)")
    parser.add_argument("--sentrix-position-column",dest="position_col",type=int,default=4,\
                        help="col in sample file where Sentrix barccode found (number cols from 0)")
    parser.add_argument("--sample-delimiter",dest="sample_delimiter",type=str,default=",",\
                        help="what separates entries in the sample file"),
    parser.add_argument("-w","--warning-only-missing-idat",dest="warning",action="store_true",default=False,\
                        help="if IDAT files are missing, report only, otherwise crash")
    args = parser.parse_args()
    return args


def parseSampleSheet(args):
    #parse the sample file to extract the IDs of the particpants and their corresponding
    #idat files. Print warning or crash if the files don't exst
    with open(args.sample) as mf:
        idats={}
        for line in mf:
            recs    = line.split(args.sample_delimiter)
            pid     = recs[args.sample_col]
            barcode = recs[args.barcode_col]
            pos     = recs[args.position_col]
            curr_fs = []
            ok= True
            warning = ""
            for colour in ["Grn","Red"]:
                f = os.path.join(args.idat,"{barcode}_{pos}_{colour}.idat".format(barcode=barcode,pos=pos,colour=colour))
                this_ok = os.access(f,os.R_OK)
                if not this_ok: warning=warning+"Warning: file {} does not exist or readable\n".format(f,EOL)
                ok = ok & this_ok
                curr_fs.append(f)
            if not ok:
                if args.warning: 
                    print(warning)
                    continue
                else:
                    sys.exit("Missing idat files: "+EOL+warning)
            idats[pid]=curr_fs
    return idats



if __name__ == '__main__':
    args = parseArguments()
    idats =  parseSampleSheet(args) 
    getIdata(args,idats)

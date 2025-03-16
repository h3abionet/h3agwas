#!/usr/bin/env python2.7

import argparse,sys
import time
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("--legend", help="legend file")
parser.add_argument("--hap", help="hap file")
parser.add_argument("--outfile", help="output file pattern (will end with .hap.gz and .legend.gz)")
parser.add_argument("--start",  type=int, help="start position")
parser.add_argument("--stop",  type=int, help="stop position")
args = parser.parse_args()

WINDOW = 500000


def extract_region(legend_file, hap_file, outfile, start, stop):
    '''
    Return: extract a region of impute reference set
    '''
    # first filter the legend fileread the legend file and negnd in the
    fin_legend = gzip.open(legend_file)
    head = fin_legend.readline()
    fin_hap = gzip.open(hap_file)

    fout_leg = gzip.open(outfile+".legend.gz","wt")
    fout_hap = gzip.open(outfile+".hap.gz","wt")
    fout_leg.write(head)

    skipped = 0
    total = 0

    for line in fin_legend:
        total += 1
        lrow = line.strip().split()
        hline = fin_hap.readline()
        if start <= int(lrow[1]) <= stop:
            fout_leg.write(line)
            fout_hap.write(hline)
        else:
            skipped+=1
        if total % 1000 == 0:
            print "\b.",

    kept = total - skipped

    fout_hap.close()
    fout_leg.close()
    print "Done. {} kept, {} skipped.".format(kept,skipped)


if __name__ == '__main__':
    extract_region(args.legend, args.hap, args.outfile, args.start-WINDOW, args.stop+WINDOW)

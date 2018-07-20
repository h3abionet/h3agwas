
from __future__ import print_function

import pandas as pd
from   Bio import SeqIO
import glob
import re
import os
import sys
import gzip
import distance
import argparse


TAB=chr(9)
EOL=chr(10)

def parseArguments():
    parser=argparse.ArgumentParser()
    parser.add_argument('reference_genome_dir', type=str, metavar='referencegenomedir'),
    parser.add_argument('strand_report', type=str, metavar='strand_report'),
    parser.add_argument('chip_manifest', type=str, metavar='chip_manifest'),
    parser.add_argument('output', type=str, metavar='output',help="output report"),
    parser.add_argument('--chrom-only',type=str,default=False,dest="chrom_only")
    parser.add_argument('--seg-len', type=int, dest='seg_len', metavar='seg_len',\
                        default = 25, help="how much matches"),
    args = parser.parse_args()
    return args


def getRef(path,build):
    match   = os.path.join(path,str(build),"genomic","*fa.gz")
    all_files = glob.glob(match)
    genome = {}
    for fname in all_files:
        m = re.search(".*.chromosome\.(\w+)\.fa.gz",fname)
        if not m:
            sys.exit("Illegal file "+fname)
        chrom = m.group(1)
        if chrom not in chroms: continue
        print(chrom,end="\t")
        seq = SeqIO.read(gzip.open(fname,"rt"),"fasta")
        genome[chrom]=seq
    print()
    return genome


def getData(directory, strand_fn, manifest_fn):
    genome = [None]*40

    for build in [36,37]:
        print("Reading in build ",build)
        genome[build]  = getRef(directory,build)

    mf     = pd.read_csv(manifest_fn,delimiter=",",dtype={"Chr":str}, skiprows=7)
    strand = \
      pd.read_csv(strand_fn,delim_whitespace=True,dtype={"Chr":str}, comment="#")
    return (genome, mf, strand)


comp = {'A':'T','C':'G','G':'C','T':'A'}


def rc(seq):
    m = ""
    for x in seq:
        m = comp.get(x,x)+m
    return m
        

def match(ref_left,probe_left,ref_right,probe_right,RC,dist):
    if RC:
        tmp = probe_left
        probe_left  = rc(probe_right)
        probe_right = rc(tmp)
    pl_len = len(probe_left)
    pl_rgt = len(probe_right)
    pleft  = probe_left[-min(seg_len,pl_len):]
    pright = probe_right[:min(seg_len,pl_rgt)]
    rleft  = ref_left[-min(seg_len,pl_len):].seq
    rright = ref_right[:min(seg_len,pl_rgt)].seq
    try:
       d1 = dist(rleft,pleft)
       d2 = dist(rright,pright)
    except ValueError:
       print("\nString length mismatch error,i,",snp["SNP_Name"])
       print("LEFT:  ",rleft,len(rleft),pleft,len(pleft),"RC=",RC)
       print("RIGHT: ",rright,len(rright),pright,len(pright))
       return 50
    return d1+d2


def  warn(warnf,chrom_num,coord,snp,direc,score):
    if score>4:
       warnf.write(TAB.join(map(str,["WARN:",chrom_num,coord,snp,direc,score]))+EOL)

def alignSNP(warnf,chromosome,coord,snp,the_snp,probe_pre,probe_pst,dist=distance.hamming):
    seq       = chromosome[coord-seg_len:coord+seg_len+1]
    base      = seq.seq[seg_len]
    ref_pre   = seq[0:seg_len]
    ref_post  = seq[seg_len+1:]
    fwd       = match(ref_pre,probe_pre,ref_post,probe_pst,False,dist)
    rev       = 2000
    score     = 1000
    align     = "-"
    if fwd<10:
       align="fwd"
       score=fwd
       warn(warnf,chrom_num,coord+1,snp["SNP_Name"],"+",fwd)
    else:
       rev=match(ref_pre,probe_pre,ref_post,probe_pst,True,dist)
       if rev<10:
          align="rev"
          score=rev
          warn(warnf,chrom_num,coord+1,snp["SNP_Name"],"-",rev)
    if the_snp in  [ "[D/I]", "[I/D]"]:
       align="fwd"
       base ="I"
    return align, base, score, fwd, rev


args = parseArguments()
if args.chrom_only:
    chroms = [ args.chrom_only ]
else:
    chroms = list(map(str, range(1,23)))+['X','Y','MT']

(genome,mf,strand) = getData(args.reference_genome_dir,args.strand_report,args.chip_manifest)




g       = open(args.output+".ref","w")
warnf   = open(args.output+".wrn","w")
errf    = open(args.output+".err","w")
seg_len = args.seg_len
for i,snp in strand.iterrows():
    align      = "fwd"
    base       = snp["Top_AlleleA"] # if can't do better
    score      = -1
    coord      = snp["Coord"]-1
    build      = int(snp['Build'])
    snp_name   = snp["SNP_Name"]
    chrom_num  = snp["Chr"]
    the_snp    = mf.loc[i]["SNP"]     
    if genome[build] == None:
        i_strand   = mf.loc[i]["IlmnStrand"]
        r_strand   = mf.loc[i]["RefStrand"]
        flip       = (i_strand == "BOT") ^ (r_strand == '-')
        if flip : base = comp.get(base,base)
        if snp_name == "rs204704":
            print(snp_name,the_snp,i_strand,r_strand,flip,base,mf.loc[i]["Name"])
        output = TAB.join(map(str, [snp_name,chrom_num,str(coord+1),base, build,"build err"]))+EOL
        errf.write(output)
        output = TAB.join(map(str, [snp["SNP_Name"],chrom_num,coord+1,base,align]))+EOL
        g.write(output)
        continue
    if chrom_num not in genome[build]:
        i_strand   = mf.loc[i]["IlmnStrand"]
        r_strand   = mf.loc[i]["RefStrand"]
        flip       = (i_strand == "BOT") ^ (r_strand == '-')
        if flip : base = comp.get(base,base)
        if snp_name == "rs204704":
            print(snp_name,the_snp,i_strand,r_strand,flip,base,mf.loc[i]["Name"])
        output = TAB.join(map(str, [snp_name,chrom_num,str(coord+1),base, build,"chrom_num_err"]))+EOL
        errf.write(output)
        output = TAB.join(map(str,[snp_name,chrom_num,coord+1,base,align]))+EOL
        g.write(output)
        continue
    top_seq = snp["Top_Seq"].upper()
    top_pre = top_seq.index("[")
    top_pst = top_seq.index("]")+1
    alleles = top_seq[top_pre+1:top_pst-1].split("/")
    probe_pre = top_seq[:top_pre]
    probe_pst = top_seq[top_pst:]
    chromosome = genome[build]["X" if chrom_num == "XY" else chrom_num]
    align, base, score, fwd, rev = alignSNP(warnf,chromosome,coord,snp,the_snp,probe_pre,probe_pst)
    if score==1000: 
        align, base, score, fwd,rev = alignSNP(warnf,chromosome,coord,snp,the_snp,probe_pre,probe_pst,distance.levenshtein)
        if score<10:
           warn(warnf,chrom_num,coord+1,snp["SNP_Name"],"BM",min(fwd,rev))
        else:
           output = TAB.join(map(str, [snp_name,chrom_num,str(coord+1),base, build,"non_align",fwd,rev]))+EOL         
           errf.write(output)
    output = TAB.join(map(str,[snp_name,chrom_num,coord+1,base,align]))+EOL
    g.write(output)
g.close()
errf.close()
warnf.close()


#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np

def get_sep(x):
    list_sep = [',', ' ', '\t']
    list_sep_num = ['COM', 'SPA', 'TAB']
    if x in list_sep:
        return x
    x = x.upper()[:3]
    if x in list_sep_num:
        return list_sep[list_sep_num.index(x)]
    print(f'\nnot found sep {x}\nexit\n')
    exit(1)

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", type=str, help="dataset file name")
parser.add_argument("-z", "--file_rsres", type=str, help="dataset file name")
parser.add_argument("-w", "--file_cross", type=str, help="dataset file name")
parser.add_argument("-o", "--out", type=str, help="output file name")
parser.add_argument("-r", "--head_rs", type=str, help="head rs")
parser.add_argument("-b", "--head_bp", type=str, help="head rs")
parser.add_argument("-c", "--head_chr", type=str, help="head rs")
parser.add_argument("-s", "--sep", type=str, help="head rs")
args = parser.parse_args()

file_i = args.file
file_rs_res = args.file_rsres
file_rs_cros = args.file_cross
file_out_rs = f"{args.out}.rs"
file_out_pos = f"{args.out}.pos"
sep = get_sep(args.sep)
chr_head = args.head_chr
bp_head = args.head_bp
rs_head = args.head_rs

data = pd.read_csv(file_i, sep=sep, header=0, quoting=3)
header_i = pd.read_csv(file_i, sep=sep, nrows=0).columns.tolist()
data.columns = header_i

def process_pos(x):
    return [min(map(int, x.split(';'))), max(map(int, x.split(';')))]

def process_chr(x):
    un = list(set(filter(lambda y: int(y) > 0, x.replace('chr', '').split(';'))))
    return un[0] if len(un) == 1 else np.nan

data['PosBeginI'], data['PosEndI'] = zip(*data[bp_head].astype(str).apply(process_pos))
data['ChrI'] = data[chr_head].astype(str).apply(process_chr)
data['Num'] = range(1, len(data) + 1)

with open(file_rs_cros, 'r') as f:
    line_cross_map = [line.strip() for line in f if len(line.split('\t')) == 5]

with open('tmplaflkjv', 'w') as f:
    f.write('\n'.join(line_cross_map))

data_cross_map = pd.read_csv('tmplaflkjv', sep='\t', header=None)
data_cross_map.columns = ["ChroNewCM", "PosDebNewCM", "PosFinNewCM", 'Strand', "Num"]

data_rs = pd.read_csv(file_rs_res, sep='\t', header=None, usecols=[0, 1, 2, 3, 4])
data_rs.columns = ["ChroNewRs", "PosNewRs", "Rs", "Ref", "Alt"]

data_m = data.merge(data_rs, left_on=rs_head, right_on="Rs", how='outer').merge(data_cross_map, left_on="Num", right_on="Num", how='outer')
data_not_found = data_m[(data_m['ChroNewCM'].isna() & data_m['ChroNewRs'].isna())]
data_not_found.to_csv(f"{args.out}.notfound.tsv", sep='\t', index=False)
data_m = data_m[~(data_m['ChroNewCM'].isna() & data_m['ChroNewRs'].isna())]
data_m.to_csv(f"{args.out}.detail.tsv", sep='\t', index=False)

data_m['ChroNew'] = data_m['ChroNewRs'].astype(str)
data_m['PosBeginNew'] = data_m['PosNewRs']
data_m['PosEndNew'] = data_m['PosNewRs']
bal = data_m['ChroNew'].isna()
data_m.loc[bal, 'ChroNew'] = data_m.loc[bal, 'ChroNewCM']
data_m.loc[bal, 'PosBeginNew'] = data_m.loc[bal, 'PosDebNewCM']
data_m.loc[bal, 'PosEndNew'] = data_m.loc[bal, 'PosFinNewCM']

data_m['ChroNew'] = data_m['ChroNew'].str.replace("chr", "")
tmp_sup1 = data_m['Num'].value_counts()
tmp_sup1 = tmp_sup1[tmp_sup1 > 1].index.tolist()
data_m_multi = data_m[data_m['Num'].isin(tmp_sup1)]
data_m = data_m[~data_m['Num'].isin(tmp_sup1)]
head_all = [col for col in data_m.columns if col not in ["Num", "ChroNewCM", "PosDebNewCM", "PosFinNewCM", "ChroNewRs", "PosNewRs"]]

data_m[head_all].to_csv(f"{args.out}.tsv", sep='\t', index=False)
data_m_multi.to_csv(f"{args.out}.multi.tsv", sep='\t', index=False)



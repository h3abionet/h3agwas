#!/usr/bin/python3
from optparse import OptionParser
import numpy as np
import pandas as pd
import gzip

parser = OptionParser()
parser.add_option("-v", "--vcf", dest="filename",
                  help="Input: Gzip compressed vcf file", type="string")
parser.add_option("-o", "--output", dest="output",
                  help="Output: csv file with SNP counts per individual", type="string")
(options, args) = parser.parse_args()

vcf_file = options.filename
output = options.output

def get_counts(vcf_file):
    total_SNPS = 0
    '''
    :param vcf_file: Whole Genome vcf file compressed with gzip
    :return: np array with SNP count, unique/private SNPs and total imputed/genotyped SNPs of each individual
    '''
    with gzip.open(vcf_file, 'rt') as fp:
        for line in fp:
            line = line.rstrip()
            if not line.startswith('#'):
                total_SNPS += 2
                columns = line.split('\t')
                num_indiv = len(columns) - 9
                indices_hetero = [i for i, item in enumerate(columns) if item.startswith('1|0') or item.startswith('0|1')]
                indices_homo = [i for i, item in enumerate(columns) if item.startswith('1|1')]
                all_indices = indices_hetero + indices_homo

                if len(all_indices) < 1:
                    continue
                elif len(indices_hetero) == 1 and len(indices_homo) == 0:
                    index = indices_hetero[0]
                    Singleton[index] += 1
                elif len(indices_homo) == 1 and len(indices_hetero) == 0:
                    index = indices_homo[0]
                    Singleton[index] += 2
                for i in indices_hetero:
                    SNP_count[i] += 1
                for j in indices_homo:
                    SNP_count[j] += 2
            elif line.startswith("##"):
                continue
            elif line.startswith('#'):
                headerlist = line.split('\t')
                SNP_count = [0] * len(headerlist)
                total = [0] * len(headerlist)
                Singleton = [0] * len(headerlist)

    for i in range(len(total)):
        total[i] = total_SNPS
    all = np.column_stack((headerlist, total, SNP_count, Singleton))
    return(all)

SNP_counts = get_counts(vcf_file)

SNP_counts_df = pd.DataFrame(SNP_counts)
SNP_counts_df = SNP_counts_df.drop(SNP_counts_df.index[0:9])
SNP_counts_df.columns = ["Indiv", "total_SNPs", "SNP_count", "Singleton_count"]

SNP_counts_df.to_csv(output, index=False, header=True)

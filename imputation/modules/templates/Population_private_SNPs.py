#!/usr/bin/python3
import numpy as np
import pandas
from optparse import OptionParser

# example for concatenating all files together into one file(INPUT):
'''
#!/usr/bin/bash
head -1 ./all.10.prep.withlabels.maf_only.with_coord > all.txt
for i in ./all.*.prep.withlabels.maf_only.with_coord
do
tail -n +2 "$i" >> all.txt
done
'''

parser = OptionParser()
parser.add_option("-i", "--input", dest="filename",
                  help="tab delimited text file with SNP MAFs of different Population", type="string")
parser.add_option("-o", "--output", dest="output",
                  help="Output: csv file with SNP counts per individual", type="string")
(options, args) = parser.parse_args()


input = options.filename
output = options.output

def get_counts(filename):
    '''
    :param filename: SNP MAFs of different Populations in a tab delimited text file:
    :return: np array with unique/private SNPs of every population
    '''

    first_line = True
    num_SNPs = 0

    with open(filename, 'rt') as fp:
        for line in fp:
            line = line.rstrip()
            list = line.split('\t')
            del list[0]
            if first_line:
                length = len(list)
                header = list
                num_unique = [0] * length
                num_snps = [0] * length
                first_line = False
            else:
                num_SNPs += 1
                index = [i for i, item in enumerate(list) if float(item) > 0]
                if len(index) > 1:
                    continue
                elif len(index) == 1:
                    num_unique[index[0]] += 1

    for i in range(len(num_snps)):
        num_snps[i] = num_SNPs

    return(np.column_stack((header, num_unique, num_snps)))


df = pandas.DataFrame(get_counts(input))
df.to_csv("Population_private_SNPs.csv", index=False, header=True)

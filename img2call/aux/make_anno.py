
import pandas as pd

import numpy as np

import sys
d = pd.read_excel("dna.xlsx")
print(d.columns)

cs = ['Row', 'Institute Sample Label', 'Sample Plate', 'Well', 'Array Info.S', 'Sentrix ID']


dtype = {'Genome' : str, 'Build' : str, 'Chr' : str }
manifest = pd.read_csv("A3.csv",skiprows=7,dtype=dtype)
freq = pd.read_csv("chip.frq",delim_whitespace=True)

mrg  = manifest.merge(freq,left_on="Name",right_on="SNP",how='outer')

mrg['isSNP'] = np.where(mrg['MAF']>0.0025,1,0)

renames = {'Chr':'chromosome', 'MapInfo':'position', 'SNP_x':'featureNames'}

mrg.rename(mapper=renames,axis=1,inplace=True)

mrg.to_csv("A3.annot",index=False)


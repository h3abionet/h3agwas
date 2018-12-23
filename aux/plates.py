# Scott Hazelhurst, University of the Witwatersrand, Johannesburg  (C) 2018 on behalf of H3ABioNet

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

sex_f = open(sys.argv[1])
the_ids=[]
for line in sex_f:
    data = line.strip().split()
    the_ids.append(data[0])



batches = pd.read_excel(sys.argv[2])
out_dir = sys.argv[3]

batches['ID'] = batches['Institute Sample Label'].apply(lambda x:x[18:])
batches['row']=batches['Well'].apply(lambda x:ord(x[0])-ord('A'))
batches['col']=batches['Well'].apply(lambda x:int(x[1:4]))
print(batches[batches['col']==0])

problems = batches[batches['ID'].isin(the_ids)]

plates=problems.groupby('Institute Plate Label')

for plate_name,plate in plates:
    plt.xlim(0,13)
    plt.ylim(-1,8)
    plt.xticks(range(13),["",1,2,3,4,5,6,7,8,9,10,11,12])
    plt.yticks(range(8),["A","B","C","D","E","F","G","H"])
    plt.scatter(plate['col'],plate['row'],color='blue')
    plt.savefig("%s/%s.png"%(out_dir,plate_name))
    plt.close()

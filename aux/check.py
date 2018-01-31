


import pandas as pd


m = pd.read_csv("/dataB/AWIGenGWAS/aux/H3Africa_2017_20021485_A3.csv",skiprows=7,usecols=["Name","IlmnStrand","RefStrand"],delimiter=",",dtype={"Chr":str})
#s = pd.read_csv("/dataB/AWIGenGWAS/aux/H3Africa_2017_20021485_A3_StrandReport_FT.txt",usecols=["SNP_Name","Forward_Allele1","Top_AlleleA"],delim_whitespace=True,comment="#",dtype={"Chr":str})

xx = (m[(m['RefStrand']=="+")&(m['IlmnStrand']=="TOP")])["Name"]
for snp in xx.values:
    print(snp)




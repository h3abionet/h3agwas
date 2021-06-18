#!/usr/bin/env python3

import sys

print('hello')

def main():

    inf = "0"
    gc10= 0
    idpat = 0
    badf = 'poorgc10.lst'

    if inf in ["","0","false","False"]:
        noGCAnalysis(badf)
        sys.exit(0)

    # if ".xls" in inf:
    #     sdf = pd.read_excel(inf)
    # else:
    #     sdf = pd.read_csv(inf,delimiter=",")

    # if "10%_GC_Score" not in sdf.columns:
    #     noGCAnalysis(outf,badf)
    #     sys.exit(0)

    # if "Institute Sample Label" not in sdf.columns:
    #     if "Sample Plate" not in sdf.columns or "Well" not in sdf.columns:
    #         sys.exit("There is no field <Institute Sample Label> in the samplesheet <%s> and I can't guess it"%inf)
    #     sdf["Institute Sample Label"] = sdf.agg('{0[Sample Plate]}_{0[Well]}'.format, axis=1)
    #     idpat="(.*)"

    # bad = sdf[sdf["10%_GC_Score"]<float(gc10)]["Institute Sample Label"].apply(extractID).values
    # num_bad = len(bad)
    # g=open(badf,"w")
    # g.writelines(map(lambda x:"%s %s"%(x,x)+EOL,bad))
    # g.close()

    # if "Call_Rate" not in sdf.columns:
    #     graphs = []
    # else:
    #     mins  = sdf[["Call_Rate","10%_GC_Score"]].min()
    #     plateg = sdf.groupby("Institute Plate Label")
    #     graphs = plotGraphs(plateg,mins)
    
    #produceTeX(outf,graphs,num_bad)


def noGCAnalysis(badf):
   g=open(badf,"w")
   g.close()

if __name__ == "__main__":
    main()
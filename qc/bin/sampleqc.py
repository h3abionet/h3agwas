#!/usr/bin/env python3

import pandas as pd
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re
import os.path

num_cols=4
num_rows=5

inf = sys.argv[1]
gc10= sys.argv[2]
idpat = sys.argv[3]
badf = sys.argv[4]
outf = sys.argv[5]

font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 12,
        }

def plotPlate(curr,plate,min_cr,min_gc):
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.rc('axes',titlesize=24)
    plt.rc('axes',labelsize=24)
    plt.scatter(plate["Call_Rate"],plate["10%_GC_Score"],color="black")
    plt.xlim(0.9*min_cr)
    plt.ylim(0.9*min_gc)
    plt.xlabel("Call Rate")
    plt.ylabel("10% GC Score")
    plt.tight_layout()
    plt.savefig(os.path.join("plates",curr))
    plt.close()
    


def plotGraphs(plateg,mins):
    names=[]
    for p_name, plate in plateg:
        curr="%s.png"%p_name
        names.append(p_name)
        plotPlate(curr,plate,*mins)
    return names

EOL=chr(10)

def header(g,graphs,num_bad):
   g.write(r"""
    \section{Plate Analysis} 

    \textbf{NB: depending on your experiment and analysis, since this
    QC step is at the plate level, this analysis may include samples
    not in your data (if the plates also include samples from other
    experiments)}.  There were %d samples with 10\%s GenCall value
    below %s on the plates -- these can be found in \url{poorgc10.lst}. Any of these that were part of the data
    are removed.  """%(num_bad,"%",gc10))
   if len(graphs)>0:
     g.write(r"""
       Analysis of call rate versus 10\%s GenCall Score
       can be found in Figure \ref{fig:ccrvgc10}. This should be used to
       set your QC parameters \url{cut_mind} and \url{gc10}.
      """ % ("%"))



def noGCAnalysis(outf,badf):
   g=open(outf,"w")
   g.write("""
    \section{Plate Analysis} 

    No plate analysis can be done because a  10\%  GenCall score can't be found in the sample sheet or 
    no sample sheet  given. 
   
   """)
   g.close()
   g=open(badf,"w")
   g.close()

       
def produceTeX(outf,graphs,num_bad):
    g = open(outf,"w")
    header(g,graphs,num_bad)

    if len(graphs)==0: 
        g.close()
        return
    g.write((r"\begin{longtable}{%s}"%(("|c"*num_cols)+r"|"))+EOL)
    g.write(r"\caption{Call rate versus 10\% GenCall score}"+EOL)
    g.write(r"\label{fig:ccrvgc10}"+EOL+r"\\\hline"+EOL)
    for r in range(0,len(graphs),num_cols):
        c_files = list(map(lambda x:r"\includegraphics[width=35mm]{plates/%s.png} "%x,graphs[r:r+num_cols]))
        c_files = c_files + ([""]*(num_cols-len(c_files)))
        p_files = list(map(lambda x:r"%s "%x,graphs[r:r+num_cols]))
        p_files = p_files + ([""]*(num_cols-len(p_files)))
        c_files = " & ".join(c_files) + r"\\"+EOL
        p_files = " & ".join(p_files) + r"\\\hline"+EOL
        g.write(c_files)
        g.write(p_files)
    g.write(r"\end{longtable}")
    g.write(r"\clearpage")
    g.close()
        
    
def extractID(x):
    m = re.search(idpat,x)
    return m.group(1)
        

if __name__ == "__main__":
    if inf in ["","0","false","False"]:
        noGCAnalysis(outf,badf)
        sys.exit(0)
    if ".xls" in inf:
        sdf = pd.read_excel(inf)
    else:
        sdf = pd.read_csv(inf,delimiter=",")
    if "10%_GC_Score" not in sdf.columns:
        noGCAnalysis(outf,badf)
        sys.exit(0)
    if "Institute Sample Label" not in sdf.columns:
        if "Sample Plate" not in sdf.columns or "Well" not in sdf.columns:
            sys.exit("There is no field <Institute Sample Label> in the samplesheet <%s> and I can't guess it"%inf)
        sdf["Institute Sample Label"] = sdf.agg('{0[Sample Plate]}_{0[Well]}'.format, axis=1)
        idpat="(.*)"
    bad = sdf[sdf["10%_GC_Score"]<float(gc10)]["Institute Sample Label"].apply(extractID).values
    num_bad = len(bad)
    g=open(badf,"w")
    g.writelines(map(lambda x:"%s %s"%(x,x)+EOL,bad))
    g.close()
    if "Call_Rate" not in sdf.columns:
        graphs = []
    else:
        mins  = sdf[["Call_Rate","10%_GC_Score"]].min()
        plateg = sdf.groupby("Institute Plate Label")
        graphs = plotGraphs(plateg,mins)
    produceTeX(outf,graphs,num_bad)
        

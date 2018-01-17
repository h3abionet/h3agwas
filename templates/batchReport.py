#!/usr/bin/env python3

from __future__ import print_function
import argparse
import sys
import re
import pandas as pd
from pandas import Series, DataFrame
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator


colours =['white','black','darkmagenta','green','cyan','maroon','green','black', 'blue', 'darkgrey', 'rosybrown', 'lightcoral', 'brown', 'firebrick', 'darkred', 'r', 'red', 'mistyrose', 'salmon', 'tomato', 'darksalmon', 'coral', 'orangered', 'lightsalmon', 'sienna', 'seashell', 'chocolate', 'saddlebrown', 'sandybrown', 'peachpuff', 'peru', 'linen', 'bisque', 'darkorange', 'burlywood',  'tan', 'papayawhip', 'moccasin', 'orange', 'wheat', 'oldlace', 'floralwhite', 'darkgoldenrod', 'goldenrod', 'cornsilk', 'gold', 'lemonchiffon', 'khaki', 'palegoldenrod', 'darkkhaki', 'ivory', 'beige', 'lightyellow', 'lightgoldenrodyellow', 'olive', 'y', 'yellow', 'olivedrab', 'yellowgreen', 'darkolivegreen', 'greenyellow', 'chartreuse', 'lawngreen', 'honeydew', 'darkseagreen', 'palegreen', 'lightgreen', 'forestgreen', 'limegreen', 'darkgreen', 'g', 'green', 'lime', 'seagreen', 'mediumseagreen', 'springgreen', 'mintcream', 'mediumspringgreen', 'mediumaquamarine', 'aquamarine', 'turquoise', 'lightseagreen', 'mediumturquoise', 'azure', 'lightcyan', 'paleturquoise', 'darkslategray', 'darkslategrey', 'teal', 'darkcyan', 'c', 'aqua', 'cyan', 'darkturquoise', 'cadetblue', 'powderblue', 'lightblue', 'deepskyblue', 'skyblue', 'lightskyblue', 'steelblue', 'aliceblue', 'dodgerblue', 'lightslategray', 'lightslategrey', 'slategray', 'slategrey', 'lightsteelblue', 'cornflowerblue', 'royalblue', 'ghostwhite', 'lavender', 'midnightblue', 'navy', 'darkblue', 'mediumblue', 'b', 'blue', 'slateblue', 'darkslateblue', 'mediumslateblue', 'mediumpurple', 'rebeccapurple', 'blueviolet', 'indigo', 'darkorchid', 'darkviolet', 'mediumorchid', 'thistle', 'plum', 'violet', 'purple', 'darkmagenta', 'm', 'fuchsia', 'magenta', 'orchid', 'mediumvioletred', 'deeppink', 'hotpink', 'lavenderblush', 'palevioletred', 'crimson', 'pink', 'lightpink']

def parseArguments():
    parser=argparse.ArgumentParser()
    parser.add_argument('base', type=str, metavar='base'),
    parser.add_argument('batch', type=str, metavar='batch'),
    parser.add_argument('batch_col', type=str, metavar='batch_col'),
    parser.add_argument('phenotype', type=str, metavar='phenotype'),
    parser.add_argument('pheno_col', type=str, metavar='pheno_col',help="phenotype colum"),
    parser.add_argument('imiss', type=str, metavar='imiss',help="indiv missingness"),
    parser.add_argument('sexcheck_report', type=str, metavar='lmiss',help="sex check report"),
    parser.add_argument('eigenvec', type=str, metavar='eigenvec',help="eigenvector"),
    parser.add_argument('genome', type=str, metavar='genome',help="genome"),
    parser.add_argument('sx_pickle', type=str, metavar='sx_pickle',help="pickle file with errors"),
    args = parser.parse_args()
    return args


TAB = chr(9)
space = chr(92)+"s+"
EOL = chr(10)
PCT = chr(37)

if len(sys.argv)<=1:
    sys.argv = ["batchReport.py","$base","$batch","$batch_col","$phenotype","$pheno_col","$imiss","$sexcheck_report","$eigenvec","$genome","$pkl"]


args = parseArguments()


blank =     """
  No batch report is done"
"""

res_text = "Table *-ref{table:batchrep:%(bname)s} on page *-pageref{table:batchrep:%(bname)s} shows the error rate by %(bname)s as found in file *-url{%(fname)s}."

res_template = """

 *-begin{table}[hb]*-centering
 *-begin{tabular}{l r r r r}*-hline
 %s & Num samples & Missing rate & %s poor & Sex Checkfail *-*-*-hline
"""



 
res_caption = """

*-caption{For each group shown, we have: the average percentage missingness in that group (100 ##*-times## the sum the number of missing calls over all individuals divided by the total number of  genotype calls over all individuals); the percentage of samples in that group with a genotyping error rate above 4 percent; and the percentage of samples in that group that fail the sex check using the standard PLINK parameters on the raw input data.}
*-label{table:batchrep:%s}

*-end{table}
"""



sub_fig_template = "*-subfloat[%s]{*-includegraphics[width=8cm]{%s}}"+EOL

rel_template = """

*-subsection{Relatedness}

Using PLINK, relatedness is computed on the using IBD and
##{*-widehat{*-pi}}## as a proxy. The data used for this analysis was
the result of the QC Phase 1 work.  All pairs of individuals with a
##{*-widehat{*-pi} *-geq ${pi_hat} }## are examined -- that individual
with the greater missingness is removed. The ##*-widehat{*-pi}## of
${pi_hat} is a parameter of the pipeline.

*-begin{itemize}
*-item %d individuals were removed because of relatedness.  The list of such individuals can be found in the file *-url{%s}.
*-end{itemize}


The overall relationship analysis is shown in the Table *-ref{tab:sex:overall} on page *-pageref{tab:sex:overall}. For each group (or if there is only one group, the group as a whole), we show the group and the number of pairs of related individuals. We also show the number of pairs with 
##*-widehat{*-pi}>0.95##, which prima facie indicates identical individuals or twins, and 
the number of pairs of individuals with ##*-widehat{*-pi}## between 0.45 and 0.95 which prima facie indicates
relationship of parent/child or siblingship.

*-begin{table}[htb]
*-begin{center}
*-begin{tabular}{l@{}Q{2cm} r r}*-hline
%s & Num *-emph{pairs} & ##*-widehat{*-pi}>0.95## &  ##*-widehat{*-pi} *-in [0.45, 0.95]## *-*-*-hline
%s*-hline
*-end{tabular}
*-end{center}
*-caption{Overall breakdown of ##*-widehat{*-pi}## anomalies: for the overall dataset and each sub-group we show the number of total number of pairs with ##*-widehat{*-pi} > $pi_hat##, the number of pairs with ##*-widehat{*-pi} > 0.95## (twins or sample duplication perhaps?), and the number of pairs with ##0.95 > *-widehat{*-pi} >0.5## (siblings or parent/child?).}
*-label{tab:sex:overall}
*-end{table}
"""

rel_mixed = "In the row, labelled *-emph{mixed}, we show the number of pairs that appear to be related to each other where the members of the pair are in different groups. Depending on sampling, related pairs coming from different groups may be unlikely and a result of sample mishandling. The pairs that are so mixed are: " 


duplicate = """
*-subsubsection*{Duplicates or very closely related individuals}

Table *-ref{table:vclose} on page *-pageref{table:vclose} shows pairs which have a ##{*-widehat{*-pi}}## value greater than 0.49 and are so very closely related. Those pairs with pi-hat values close to 1 are likely to be duplicates or identical twins (or a sample mishandling error). The groups/phenotypes/sites of the individuals are shown in parentheses.
"""

# no longer used --kept in case we decide to add back
duplicate_table  = """
*-begin{longtable}{l l r}*-hline
*-caption{Pairs of individuals with very high PI-hat values. These should be investigated}
*-label{table:vclose}
*-endhead
Individual 1 & Individual 2 & ##{*-widehat{*-pi}}## *-*-*-hline
%s
*-end{longtable}

"""


def getHeader():
  return """
*-section{Batch report}
*-label{sec:batch}

The batch quality results shown here (missingness and sex-check failures) is based on the raw input data (other than duplicate SNPs being deleted) as this gives insight into the quality of the raw data differences between batches and sub-batches. The PC analysis is based on the refined data.

*-subsection{Overall missingness and sex anomaly report}

"""


def poorFn(x):
    return 100*len(list(filter(lambda f: f>0.04,x)))/len(x)

def sexCheckProblem(x):
    return 100*len(list(filter(lambda p:p=='PROBLEM',x)))/len(x)

def getAnomSex(pheno_col,gname,gdf):
    anom = gdf[gdf['STATUS']=='PROBLEM']
    fname = "anom-{}-{}.lst".format(pheno_col,gname)
    f = open(fname,"w")
    for i,r in gdf.iterrows():
        f.write(" ".join(i)+EOL)
    f.close()
    return fname

def miss_vals(ifrm,pfrm,pheno_col,sexcheck_report):
    def group_fn(x):
        return pfrm.ix[x][pheno_col]
    g  = ifrm.groupby(group_fn)
    num_samples = g['N_MISS'].count()
    ave_miss    = 100*g['N_MISS'].sum()/g['N_GENO'].sum()
    num_poor_i  = g[['F_MISS']].agg(poorFn)
    sxfrm       = pd.read_csv(sexcheck_report,delim_whitespace=True,index_col=[0,1])
    g = sxfrm.groupby(group_fn)
    problems = g[['STATUS']].agg(sexCheckProblem)
    sex_report=EOL+EOL
    #"Any samples with anomalous sex status can be found in the following files: "
    #glist = []
    #for gname, gdf in g:
    #    glist.append(getAnomSex(pheno_col,gname,gdf))
    #sex_report = sex_report + ", ".join(glist)+"."+EOL+EOL
    return (num_samples,ave_miss,num_poor_i,problems,sex_report)



def showResult(colname,num_samples,ave_miss,num_poor_i,problems,sex_report):
   t = res_template%(colname,"*-"+PCT)
   for r in ave_miss.index:
      res = ( r, num_samples.ix[r],ave_miss.ix[r],num_poor_i.ix[r]['F_MISS'],problems.ix[r]['STATUS'])
      t = t + "%s & %d & %5.2f & %5.2f & %5.2f *-*-"%res + EOL
   t=t+"*-hline"+EOL+"*-end{tabular}"+(res_caption%colname)
   t=t+sex_report
   return t.replace("*-",chr(92))






def plotPCs(base,eigs,pfrm,pheno_col,batch,batch_col):
    def group_fn(x):
       try:
          return pfrm.ix[x][args.pheno_col]
       except KeyError:
          return 0
    g = eigs.groupby(group_fn)
    matplotlib.rcParams.update({'font.size': 12})
    allplots=[]
    for site, df in g:
        fig = plt.figure(figsize=(6,6))
        fig, ax = plt.subplots()
        locator = MaxNLocator(nbins=4) 
        ax.xaxis.set_major_locator(locator)
        ax.scatter(df['PC1'],df['PC2'],s=1,c=list(map(lambda x:colours[x],batch[batch_col])),label=[1,2,3])
        ax.legend(scatterpoints=1)
        recs=[]
        classes=[]
        for i in batch[batch_col].unique():
            recs.append(mpatches.Rectangle((0,0),1,1,fc=colours[i]))
            classes.append(batch_col+str(i))
        plt.legend(recs,classes,loc=4)
        plt.xlabel("PC1")
        plt.ylabel("PC2")
        plt.tight_layout()
        fn = "%s-%s.pdf"%(base,site)
        allplots.append((fn, pheno_col+" "+str(site)))
        plt.savefig(fn,type="pdf")
    return allplots

       
fig_template= """
*-begin{figure}[ht]
%s
*-caption{Principal component analysis, using  the batch of the samples as the label.}
*-label{fig:batch:pca%s}

*-end{figure}
"""



def showFigs(figs):
    num_cols=2
    num_sub_figs_per_fig=6
    if len(figs)==1:
       inner="*-includegraphics[width=8cm]{%s}"%figs[0][0]
    else:
       inner=template=""
       caps = list('abcdefghijklmnopqrtuvwxyzABCD')
       curr_cap=0
       for (i,(fig,label)) in enumerate(figs):
          inner=inner+sub_fig_template%(label,fig)
          if i%num_cols==1:
              inner=inner+"*-*-"
          if (i+1)%num_sub_figs_per_fig==0:
              template= template + fig_template%(inner,caps[curr_cap])
              inner=""
              curr_cap=curr_cap+1
       if (i+1)%num_sub_figs_per_fig !=0:
           template= template + fig_template%(inner,caps[curr_cap])
       else:
           curr_cap=curr_cap-1
       if curr_cap==0:
          ref="Figure *-ref{fig:batch:pcaa} shows a principal component analysis by \
               batch and phenotype."
       elif curr_cap==1:
          ref="*-noindent Figures *-ref{fig:batch:pcaa} and *-ref{fig:batch:pcab} show a principal component analysis by batch and phenotype."
       else:
          ref="*-noindent Figures *-ref{fig:batch:pcaa}--*-ref{fig:batch:pca%s} show a principal component analysis by batch and phenotype."%(caps[curr_cap])

       ref = "*-subsection{Principal component anlysis}"+EOL+EOL+\
             ref+""" In the PC analysis, the principal components were computed from all
             samples in the data together, and in each figure we extract the
             samples for that analysis. The labels and scale of the axes may
             differ, but the coordinates are the same: the position (0,0) is the
             same in all sub-figures."""
    return (EOL+ref+EOL+template)



def pstr(p):
    return str(p[0])+":"+str(p[1])

def getVClose(gfrm,pfrm,pheno_col):
    vclose = gfrm[gfrm['PI_HAT']>0.49][["FID1","IID1","FID2","IID2","PI_HAT"]]
    curr = TAB.join(["FID1","IID1",pheno_col+"1","FID2","IID2",pheno_col+"2","PI_HAT"])+EOL
    for i,row in vclose.iterrows():
        curr = curr+TAB.join([row[0],row[1],pfrm.loc[(row[0],row[1])][pheno_col],\
                              row[2],row[3],pfrm.loc[(row[2],row[3])][pheno_col],\
                              str(row[4])])+EOL
    return curr

def getRelatedPairs(pfrm,pheno_col,genome):
    our_types = { 'PI_HAT' : float }
    gfrm = pd.read_csv(genome,delim_whitespace=True, dtype=our_types)
    group = {}
    mixed = {}
    ident = {}
    sib   = {}
    num_mixed = ident_mixed = sib_mixed = 0
    for i,pair in gfrm.iterrows():
        id1 = pair[["FID1","IID1"]].values
        id2 = pair[["FID2","IID2"]].values
        this_g1 = pfrm[pheno_col].loc[tuple(id1)]
        this_g2 = pfrm[pheno_col].loc[tuple(id2)]
        group[" ALL"] = group.get(" ALL",0)+1
        if this_g1 == this_g2:
            group[this_g1] = group.get(this_g1,0)+1
            if pair['PI_HAT']>0.95:
                ident[this_g1]=ident.get(this_g1,0)+1
                ident[" ALL"]=ident.get(" ALL",0)+1                
            elif pair['PI_HAT']>0.45:
                sib[this_g1]=sib.get(this_g1,0)+1
                sib[" ALL"]=sib.get(" ALL",0)+1                
        else:
            if pair['PI_HAT']>0.9:
                ident_mixed = ident_mixed+1
                ident[" ALL"]=ident.get(" ALL",0)+1                
            elif pair['PI_HAT']>0.45:
                sib_mixed   = sib_mixed+1
                sib[" ALL"]=sib.get(" ALL",0)+1                
            pair = tuple(sorted((this_g1,this_g2)))
            if pair in mixed:
                mixed[pair].append((id1,id2))
            else:
                mixed[pair]=[(id1,id2)]
            num_mixed = num_mixed+1
    rows = ""

    keys = sorted(group.keys())
    for k in keys:
        rows = rows + "%s & %d & %d & %d *-*-"%(k,group[k],ident.get(k,0),sib.get(k,0))+EOL
    rel_text=""
    if num_mixed > 0:
        rows = rows + "mixed & %d & %d & %d*-*-"%(num_mixed,ident_mixed,sib_mixed)+EOL
        if num_mixed < 10:
            rel_text = "{*-footnotesize"
            for k, prs in mixed.items():
                rel_text = rel_text + EOL + EOL + (" Mixed group (%s, %s): "%k)
                for (p1, p2) in prs:
                    rel_text=rel_text+" "+pstr(p1)+" "+pstr(p2)+"; "
            rel_text=rel_mixed+rel_text+"}"+EOL
    num_rel = len(open("$rel_indivs").readlines())
    rel_file="%s-reltable.csv"%(args.base)
    text=rel_template%(num_rel,rel_file,pheno_col,rows)+rel_text
    vclose  = getVClose(gfrm,pfrm,pheno_col)
    g=open(rel_file, "w")
    g.write(vclose)
    g.close()
    return text



    

det_sex_analysis = """
*-subsection{Detailed Sex Check Analysis}

The this section we show detailed analysis of sex check errors. The purpose of this analysis is to help
identify trends between sub-groups, as well as possible labelling and sampling errors. In this analysis,
we use PLINK to analyse the non-recombining regions of the X-chromosome, and in particular its computation 
of the inbreeding co-efficient of the X-chromosome. If the F statistic is greater than $f_low_male, PLINK
infers that the sample is male; if it is less than $f_hi_female, it infers that the sample is female.

There are two types of errors that can happen. *-emph{Soft} errors are
those cases where an individual is slightly above or below the stated
threshholds. These may not be errors -- since the F cut-off values are
arbitrary, a too strict F-value may be chosen, or there may be unusual
patterns within the individuals studied. How these samples should be
treated will require some thought, but these are not a sign of
problems of the experimental protocol, and not by itself probably a sign of problems with
genotyping or DNA quality errors. We define a soft error as indvidual having an ##F## value beteen 0.34 and 0.66.

*-emph{Hard} errors are cases where the sex in the manifest/fam file
is markedly different from what the F-statistic predicts (e.g., the
F-statistic says 0.996 and we have this as a female). While there may
be very unusual cases were this occurs, many such examples in the data
are likely to be a sign of sample handling errors. Great care needs to
be taken with this.

A summary of the detailed analysis is shown below. In the output, the
file *-url{%s} is a CSV file that has sample ID, per-individual
missingness rate in the raw data, the %s status, and whether that
sample is hard (H) or soft error (S) for given tolerated per-individual
genotyping error rates on the X-chromosome. A `-' indicates that the indivdual is filtered out
at that rate of missingness This can be used to assess
whether anomalous sex results are due to poor genotyping rates in
individuals.


"""


det_table="""*-begin{table}[htb]*-centering

*-begin{tabular}{%s}*-hline
%s
%s*-*-*-hline
%s*-hline
*-end{tabular}


*-caption{Summary of anomalous sex calls. The results are shown for
different values *-emph{mind} the maximum per-sample error rate
tolerated in the X-chromosome (PLINK parameter). *-emph{mind}==1 means
all SNPs included. *-emph{Tot} shows the number of individuals
included with the given *-emph{mind} value. *-emph{HErr} is the number of hard errors. 
*-emph{SErr} is the number of soft errors.}

*-label{tab:sxdet:%s} *-end{table}

"""

def detSexHeader(sxAnalysis,pheno_col):
    num_bands = len(sxAnalysis.columns) # missing rates supported
    tab_spec= "l"+(" r@{*-phantom{.}}r@{*-phantom{.}}r"*num_bands)
    header1 = pheno_col 
    for n in sxAnalysis.columns:
       header1=header1+" & *-multicolumn{3}{c}{mind=%s}"%str(n)
    header1 = header1+"*-*-"
    header2 = "& Tot & SErr & HErr"*num_bands+"*-*-*-hline"
    return (tab_spec,header1,header2)

# provide the details of a specific group
def detSexGroup(sfrm,gname):
    tot   = sfrm.count()-sfrm.isnull().sum()
    herrs = (sfrm=='H').sum()
    serrs = (sfrm=='S').sum()
    line = "%8s"%gname 
    for rate in sfrm.columns:
        line = line+" & %d & %d & %d "%(tot[rate],serrs[rate],herrs[rate])
    line=line+"*-*-"+EOL
    return line
   

def xstr(m):
    if type(m) == type(0.5) and np.isnan(m):
        return "-"
    else:
        return str(m)

def dumpMissingSexTable(fname, ifrm,sxAnalysis,pfrm):
    g=open(fname,"w")
    g.write(TAB.join(["FID","IID"])+TAB+TAB.join(map(str,sxAnalysis.columns))+EOL)
    for i, row in ifrm.iterrows():
        output = TAB.join(map(str, [*i,"%5.3f"%row['F_MISS'],pfrm.loc[i]]))+\
                 TAB+TAB.join(map(xstr,sxAnalysis.loc[i]))+EOL
        g.write(output)
    g.close()



def detailedSexAnalysis(pfrm,ifrm,sxAnalysisPkl,pheno_col):
    def group_fn(x):
        return pfrm.ix[x][pheno_col]
    sex_fname = args.base+"_missing_and_sexcheck.csv"
    sxAnalysis = pd.read_pickle(sxAnalysisPkl)
    dumpMissingSexTable(sex_fname,ifrm,sxAnalysis,pfrm[pheno_col])
    header = detSexHeader(sxAnalysis,pheno_col)
    tbl    = detSexGroup(sxAnalysis,"overall")
    g = sxAnalysis.groupby(group_fn)
    for grpname, gg in g:
        tbl = tbl+detSexGroup(gg,grpname)

    return \
        det_sex_analysis%(sex_fname,pheno_col) +\
        det_table%(header+(tbl,pheno_col))

       
def backslashify(text):
    return text.replace("*-",chr(92)).replace("##",chr(36))


text = getHeader()+EOL+EOL


#if not(args.batch == "0" and args.phenotype == "0"):
ifrm = pd.read_csv(args.imiss,delim_whitespace=True,index_col=[0,1])

if args.batch != "0":
    text = text+res_text%({'bname':args.batch_col,'fname':args.batch})
    bfrm = pd.read_csv(args.batch,delim_whitespace=True,index_col=[0,1])
    problems = ifrm.index.difference(bfrm.index).to_series().values
    if len(problems)>0:
        print("Problem there IDs in the genotype data that are not in the batch data")
        print("You need to fix one way or the other")
        print("The IDs are")
        for p in problems:
            print(p)
        sys.exit(-1)
    result = miss_vals(ifrm,bfrm,args.batch_col,args.sexcheck_report)
    text = text+showResult(args.batch_col,*result)
if args.phenotype !="0":
    pfrm = pd.read_csv(args.phenotype,delim_whitespace=True,index_col=[0,1])
    text = text+res_text%({'bname':args.pheno_col,'fname':args.phenotype})
    result = miss_vals(ifrm,pfrm,args.pheno_col,args.sexcheck_report)
    text = text+showResult(args.pheno_col,*result)

col_names=['FID','IID']+list(map(lambda x: "PC%d"%x,range(1,21)))

eigs=pd.read_csv(args.eigenvec,delim_whitespace=True,header=None,names=col_names,index_col=[0,1])

if args.batch in ["0",False,"False","FALSE","nil","false",0,None,""]::
    args.batch_col = 'batch'
    bfrm = DataFrame([1]*len(eigs),index=eigs.index,columns=['batch'])
if args.phenotype in ["0",False,"False","FALSE","nil","false",0,None,""]:
    args.pheno_col = 'all'
    pfrm = DataFrame([1]*len(eigs),index=eigs.index,columns=['all'])

text = text + detailedSexAnalysis(pfrm,ifrm,args.sx_pickle,args.pheno_col)


m = re.search("(.*).eigenvec",args.eigenvec)
base = m.group(1)
figs = plotPCs(base,eigs,pfrm,args.pheno_col,bfrm,args.batch_col)
figs_text = showFigs(figs)
related_text = getRelatedPairs(pfrm,args.pheno_col,args.genome)
text = backslashify(text+related_text+figs_text)
g=open("%s-batch.tex"%args.base,"w")
g.write(text)
g.close()
   

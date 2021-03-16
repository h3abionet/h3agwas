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

def debug(*args):
    return
    print(*args)


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

null_file = "emptyZ0"

if len(sys.argv)<=1:
    sys.argv = ["batchReport.py","$base","$batch","$batch_col","$phenotype","$pheno_col","$imiss","$sexcheck_report","$eigenvec","$genome","$pkl"]


args = parseArguments()


idtypes = dict(map(lambda x: (x,str),\
                   ["FID","IID","FID1","IID1","FID2","IID2",str(args.pheno_col), str(args.batch_col)]))
idtypes['PI_HAT']=float


blank =     """
  No batch report is done"
"""

def getCsvI(fn,names=None):
    ''' read a CSV file with IDs as index '''
    frm = pd.read_csv(fn,delim_whitespace=True,names=names,dtype=idtypes)
    frm = frm.set_index(['FID','IID'])
    # nb bug/feature in pandas, if a column is an index then dtype is not applied to it
    # can have numbers as IDs but want them as string
    return frm


res_template = """

 *-begin{table}[hb]*-centering
 *-begin{tabular}{l r r r r}*-hline
 *-url{%s} & Num samples & Missing rate & *-url{%s} poor & Sex Checkfail *-*-*-hline
"""



 
res_caption = """

*-caption{For each group shown, we have: the average percentage missingness in that group (100 ##*-times## the sum the number of missing calls over all individuals divided by the total number of  genotype calls over all individuals); the percentage of samples in that group with a genotyping error rate above 4 percent; and the percentage of samples in that group that fail the sex check using the standard PLINK parameters on the raw input data.}
*-label{table:batchrep:%s}

*-end{table}
"""



sub_fig_template = "*-subfloat[%s]{*-includegraphics[width=8cm]{%s}}"+EOL

none_rel_template = """

*-clearpage

*-subsection{Relatedness}

There are no related pairs of individuals (##*-widehat{*-pi} > ${pi_hat}##, the parameter of the pipeline).


"""

rel_template = """

*-clearpage

*-subsection{Relatedness}

Using PLINK, relatedness is computed on using IBD with
##{*-widehat{*-pi}}## as a proxy. The data used for this analysis was
the result of the QC Phase 1 work.  The ##*-widehat{*-pi}## of
${pi_hat} is a parameter of the pipeline.

All pairs of individuals with a
##{*-widehat{*-pi} *-geq ${pi_hat} }## are examined -- we try to remove as few individuals as possible by
removing first those who are related to multiple people (e.g. if A is a cousin of B and C, B and C may not be related so it makes sense to remove A rather than B and C).

*-begin{itemize}
*-item There were %(numpairs)d pairs of individuals over the cut-off. These can be found in the PLINK report *-url{%(all_pi_hat)s}.
*-item %(num_rem)d individuals were removed because of relatedness.  The list of such individuals can be found in the file *-url{${rem_indivs}}. (One individual in each pair is removed;  any individuals in a pair with relatedness strictly greater than $super_pi_hat are removed.)
*-item The individuals with ##{*-widehat{*-pi} *-geq 0.45 }## can be found in *-url{%(vclose)s}.
*-end{itemize}


The overall relationship analysis is shown in the Table *-ref{tab:sex:overall} on page *-pageref{tab:sex:overall}. For each group (or if there is only one group, the group as a whole), we show the group and the number of pairs of related individuals. We also show the number of pairs with 
##*-widehat{*-pi}>0.95##, which prima facie indicates identical individuals or twins, and 
the number of pairs of individuals with ##*-widehat{*-pi}## between 0.45 and 0.95 which prima facie indicates relationship of parent/child or siblingship.

*-begin{table}[htb]
*-begin{center}
*-begin{tabular}{l@{}Q{2cm} r r r}*-hline
*-url{%(pheno_col)s} & Num *-emph{pairs} & ##*-widehat{*-pi}>0.95## &  ##*-widehat{*-pi} *-in [0.45, 0.95]## & ##<0.45## *-*-*-hline
%(rows)s*-hline
*-end{tabular}
*-end{center}
*-caption{Overall breakdown of ##*-widehat{*-pi}## anomalies: for the overall dataset and each sub-group we show the number of total number of pairs with ##*-widehat{*-pi} > $pi_hat##, the number of pairs with ##*-widehat{*-pi} > 0.95## (twins or sample duplication perhaps?), and the number of pairs with ##0.95 > *-widehat{*-pi} >0.45## (siblings or parent/child?).}
*-label{tab:sex:overall}
*-end{table}
"""

rel_mixed = "In the row, labelled *-emph{mixed}, we show the number of pairs that appear to be related to each other where the members of the pair are in different groups. Depending on sampling, related pairs coming from different groups may be unlikely and a result of sample mishandling. The pairs that are so mixed are: " 


duplicate = """
*-subsubsection*{Duplicates or very closely related individuals}

Table *-ref{table:vclose} on page *-pageref{table:vclose} shows pairs which have a ##{*-widehat{*-pi}}## value greater than 0.45 and are so very closely related. Those pairs with pi-hat values close to 1 are likely to be duplicates or identical twins (or a sample mishandling error). The groups/phenotypes/sites of the individuals are shown in parentheses.
"""

# no longer used --kept in case we decide to add back
duplicate_table  = """
*-begin{longtable}{l l r}*-hline
*-caption{Pairs of individuals with very high PI-hat values. These should be investigated}
*-label{table:vclose}
*-endhead
Individual 1 & Individual 2 & ##{*-widehat{*-pi}}## *-*-*-hline
*-url{%s}
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
    g  = pd.merge(pfrm,ifrm,left_index=True,right_index=True,how='inner').groupby(pheno_col)
    num_samples = g['N_MISS'].count()
    ave_miss    = 100*g['N_MISS'].sum()/g['N_GENO'].sum()
    num_poor_i  = g[['F_MISS']].agg(poorFn)
    if "$extrasexinfo" == "--must-have-sex":
        sxfrm       = getCsvI(sexcheck_report)
        g = pd.merge(pfrm,sxfrm,left_index=True,right_index=True,how='inner').groupby(pheno_col)
        problems = g[['STATUS']].agg(sexCheckProblem)
    else:
        problems = None
    sex_report=EOL+EOL
    #"Any samples with anomalous sex status can be found in the following files: "
    #glist = []
    #for gname, gdf in g:
    #    glist.append(getAnomSex(pheno_col,gname,gdf))
    #sex_report = sex_report + ", ".join(glist)+"."+EOL+EOL
    return (num_samples,ave_miss,num_poor_i,problems,sex_report)



def showResult(colname,num_samples,ave_miss,num_poor_i,problems,sex_report):
   t = res_template%(colname,"*-"+PCT)
   if type(problems) == type(None): 
       template="*-url{%s} & %d & %5.2f & %5.2f & n/a *-*-"
   else:
       template="*-url{%s} & %d & %5.2f & %5.2f & %5.2f *-*-"
   for r in ave_miss.index:
      if type(problems) == type(None):
         res = ( r, num_samples.loc[r],ave_miss.loc[r],num_poor_i.loc[r]['F_MISS'])
      else:
         res = ( r, num_samples.loc[r],ave_miss.loc[r],num_poor_i.loc[r]['F_MISS'],problems.loc[r]['STATUS'])
      t = t + template%res + EOL
   t=t+"*-hline"+EOL+"*-end{tabular}"+(res_caption%colname)
   t=t+sex_report
   return t.replace("*-",chr(92))






def plotPCs(base,eigs,pfrm,pheno_col,batch,batch_col):
    debug("In plotPCs")
    debug("eigs",eigs.index.names," columns", eigs.columns)  #DBG
    debug("pfrm",pfrm.index.names," columns",pfrm.columns)  #DBG
    debug("batch",batch.index.names," batch",batch.columns) #DBG
    all_f = eigs.join(pfrm)
    g = all_f.groupby(args.pheno_col)
    matplotlib.rcParams.update({'font.size': 12})
    batch_values = batch[batch_col].unique()
    our_colours  = dict(zip(batch_values,colours[1:1+len(batch_values)]))
    allplots=[]
    batch['colour']=batch.apply(lambda x:our_colours[x[batch_col]],axis=1)
    for site, df in g:
        debug("plotPCs loop",site,df.head()) #DBG
        df=df.merge(batch,left_index=True,right_index=True,how='inner')
        fig = plt.figure(figsize=(6,6))
        fig, ax = plt.subplots()
        locator = MaxNLocator(nbins=4) 
        ax.xaxis.set_major_locator(locator)
        ax.scatter(df['PC1'],df['PC2'],s=1,\
                   c=df['colour'],\
                   label=[1,2,3])
        ax.legend(scatterpoints=1)
        recs=[]
        classes=[]
        for i in batch[batch_col].unique():
            recs.append(mpatches.Rectangle((0,0),1,1,fc=our_colours[i]))
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
    ref = template = ""
    if len(figs)==1:
       template=fig_template%("*-includegraphics[width=8cm]{%s}"%figs[0][0],"onlypca")+EOL+EOL
       ref="""In the PC analysis, the principal components were computed from all
             samples in the data together, and in the figure  we extract the
             samples for that analysis. See Figure *-ref{fig:batch:pcaonlypca}."""
    else:
       inner=template=""
       caps = list('abcdefghijklmnopqrtuvwxyzABCD')
       curr_cap=0
       for (i,(fig,label)) in enumerate(figs):
          inner=inner+sub_fig_template%(label.replace('_',' '),fig)
          if i%num_cols==1:
              inner=inner+"*-*-"
          if (i+1)%num_sub_figs_per_fig==0:
              template= template + fig_template%(inner.replace('_',' '),caps[curr_cap])
              inner=""
              curr_cap=curr_cap+1
       if (i+1)%num_sub_figs_per_fig !=0:
           template= template + fig_template%(inner.replace('_',' '),caps[curr_cap])
       else:
           curr_cap=curr_cap-1
       if curr_cap==0:
          ref="Figure *-ref{fig:batch:pcaa} shows a principal component analysis by \
               batch and phenotype."
       elif curr_cap==1:
          ref="*-noindent Figures *-ref{fig:batch:pcaa} and *-ref{fig:batch:pcab} show a principal component analysis by batch and phenotype."
       else:
          ref="*-noindent Figures *-ref{fig:batch:pcaa}--*-ref{fig:batch:pca%s} show a principal component analysis by batch and phenotype."%(caps[curr_cap])

       ref = """ In the PC analysis, the principal components were computed from all
             samples in the data together, and in the figure(s) we extract the
             samples for that analysis. The labels and scale of the axes may
             differ, but the coordinates are the same: the position (0,0) is the
             same in all sub-figures."""
    ref="*-subsection{Principal component anlysis}"+EOL+EOL+ref
    return (EOL+ref+EOL+template)



def pstr(p):
    return str(p[0])+":"+str(p[1])

def getVClose(gfrm,pfrm,pheno_col):
    vclose = gfrm[gfrm['PI_HAT']>0.45][["FID1","IID1","FID2","IID2","PI_HAT"]]
    curr = TAB.join(["FID1","IID1",pheno_col+"1","FID2","IID2",pheno_col+"2","PI_HAT"])+EOL
    for i,row in vclose.iterrows():
        #if re.search("_replicate_",row[0]+row[2]): continue
        #if re.search("_MSK_",row[1]+row[3]): continue
        if row[1][-4:-1]=="DUP" or row[3][-4:-1]=="DUP": continue
        pairA=pfrm.loc[row[0],row[1]][pheno_col]
        pairB=pfrm.loc[row[2],row[3]][pheno_col]
        if type(pairA)!=str: # hack to handle fake phenos -- this should be fixed upstream
            pairA = pairA.values[0]
            pairB = pairB.values[0]
        curr = curr+TAB.join([row[0],row[1],str(pairA),\
                              row[2],row[3],str(pairB),\
                              str(row[4])])+EOL
    return curr


def getGroupNum(pfrm,pheno_col,fid,iid):
   this_g = pfrm[pheno_col].loc[fid,iid]
   if 'values' in dir(this_g): this_g=this_g.values[0]
   return this_g

def getRelatedPairs(pfrm,pheno_col,genome):
    gfrm = pd.read_csv(genome,delim_whitespace=True, dtype=idtypes)
    group = {}
    mixed = {}
    ident = {}
    sib   = {}
    num_mixed = ident_mixed = sib_mixed = 0
    for i,pair in gfrm.iterrows():
        id1 = pair[["FID1","IID1"]].values
        id2 = pair[["FID2","IID2"]].values
        if "DUP" == id2[1][-4:-1] or "DUP"==id1[1][-4:-1]:continue
        #if re.search("_replicate_",id1[0]+id2[0]): continue
        #if re.search("_MSK_",id1[1]+id2[1]): continue
        this_g1 = getGroupNum(pfrm,pheno_col,id1[0],id1[1])
        this_g2 = getGroupNum(pfrm,pheno_col,id2[0],id2[1])
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
    if len(keys) == 0:
        return none_rel_template
    if len(keys)==2:
        keys=[" ALL"]
    for k in keys:
        rest = group[k]-ident.get(k,0)-sib.get(k,0)
        rows = rows + "*-url{%s} & %d & %d & %d & %d*-*-"%(k,group[k],ident.get(k,0),sib.get(k,0),rest)+EOL
    rel_text=""
    if num_mixed > 0:
        rest = num_mixed-ident_mixed-sib_mixed
        rows = rows + "mixed & %d & %d & %d & %d*-*-"%(num_mixed,ident_mixed,sib_mixed,rest)+EOL
        if num_mixed < 10:
            rel_text = "{*-footnotesize"
            for k, prs in mixed.items():
                rel_text = rel_text + EOL + EOL + (" Mixed group (%s, %s): "%k)
                for (p1, p2) in prs:
                    rel_text=rel_text+" "+pstr(p1)+" "+pstr(p2)+"; "
            rel_text=rel_mixed+rel_text+"}"+EOL
    num_rem = len(open("$rem_indivs").readlines())
    rel_file="%s-reltable.csv"%(args.base)
    rdict = { 'numpairs' : group[" ALL"], 'num_rem':num_rem, 'all_pi_hat':genome, 'vclose':rel_file, \
              'pheno_col':pheno_col, 'rows':rows }
    text=rel_template%rdict+rel_text
    # (group[" ALL"],num_rem,grel_file,pheno_col,rows)
    vclose  = getVClose(gfrm,pfrm,pheno_col)
    g=open(rel_file, "w")
    g.write(vclose)
    g.close()
    return text



    

det_sex_analysis = """
*-subsection{Detailed Sex Check Analysis}

This section shows a detailed analysis of sex check anomalies and/or
unusual patterns in the X-chromosome. The purpose of this analysis is to help identify
trends between sub-groups, as well as possible labelling and sample handling
errors. The term *-emph{anomaly} or *-emph{error} is used to label individuals where the
sex of the individual as described in the manifest does not *-emph{stricly} match analysis of the X-chromosome and
so for QC should be considerd further.


In this analysis, we use PLINK to analyse the non-recombining
regions of the X-chromosome, and in particular its computation of the
inbreeding co-efficient of the X-chromosome. If the ##F## statistic is
greater than $f_lo_male, PLINK infers that the sample is male; if it
is less than $f_hi_female, it infers that the sample is female.

Reminder: the checking of sex on the raw data before any other QC is shown in Section 1. The rest of this section analyses the data after basic QC on genotype has been done and so differences may be seen.

There are two types of apparent anomaly that can happen. *-emph{Soft} anomalies are
those cases where an individual is slightly above or below the stated
threshholds. These may not be sample handling errors -- since the ##F## cut-off values are
arbitrary, a too strict ##F##-value may be chosen, or there may be uncommon
patterns within the individuals studied. How these samples should be
treated will require some thought, but these are not a sign of
problems of the experimental protocol, and not by itself probably a sign of problems with
genotyping or DNA quality errors. We define a soft anomaly as an indvidual having an ##F## value between 
$f_hi_female and $f_lo_male, which is possible but unusual.


*-emph{Hard} errors are cases where the sex in the manifest/fam file
is markedly different from what the F-statistic predicts (e.g., the
F-statistic says 0.996 and we have this as a female). While there may
be very unusual cases were this occurs, many such examples in the data
are likely to be a sign of sample handling errors. Great care needs to
be taken with this.

A summary of the detailed analysis is shown below. In the output, the
file *-url{%s} is a CSV file that has sample ID, per-individual
missingness rate in the raw data (*-verb!F_MISS! is the individual
missingness rate *-emph{not} the F-statistic), the *-url{%s} status, and
whether that sample is hard (H) or soft error (S) for given tolerated
per-individual genotyping error rates on the X-chromosome. A `-'
indicates that the indivdual is filtered out at that rate of
missingness This can be used to assess whether anomalous sex results
are due to poor genotyping rates in individuals.


"""


det_table="""*-begin{table}[htb]*-centering

*-begin{tabular}{%s}*-hline
%s
%s*-*-*-hline
%s*-hline
*-end{tabular}


*-caption{ of anomalous sex calls. The results are shown for
different values *-emph{mind} the maximum per-sample error rate
tolerated in the X-chromosome (PLINK parameter). *-emph{mind}==1 means
all SNPs included. *-emph{Tot} shows the number of individuals
included with the given *-emph{mind} value. *-emph{HErr} is the number of hard errors. 
*-emph{SAnm} is the number of soft apparent anomalies or unusual results.}

*-label{tab:sxdet:%s} *-end{table}

"""

def detSexHeader(pfrm,missing_sex_columns):
    num_bands = len(missing_sex_columns) # missing rates supported
    tab_spec= "l"+(" r@{*-phantom{.}}r@{*-phantom{.}}r"*num_bands)
    header1 = args.pheno_col.replace('_',' ')
    for n in missing_sex_columns:
       header1=header1+" & *-multicolumn{3}{c}{mind=%s}"%str(n)
    header1 = header1+"*-*-"
    header2 = "& Tot & SAnm & HErr"*num_bands+"*-*-*-hline"
    return (tab_spec,header1,header2)

# provide the details of a specific group
def detSexGroup(sfrm,gname,missing_sex_columns):
    line = "%8s"%gname 
    for rate in missing_sex_columns:
        tot   = sfrm[rate].count()-sfrm[rate].isnull().sum()
        herrs = (sfrm[rate]=='H').sum()
        serrs = (sfrm[rate]=='S').sum()
        line = line+" & %d & %d & %d "%(tot,serrs,herrs)
    line=line+"*-*-"+EOL
    return line
   

def xstr(m):
    if type(m) == type(0.5) and np.isnan(m):
        return "-"
    else:
        return str(m)

def dumpMissingSexTable(fname, pfrm, sex_missing_cols):
    if args.batch_col+"-b" in pfrm.columns.values: # could have same col in batch and pheno file
        bcol = args.batch_col+"-b"
    else:
        bcol = args.batch_col
    cols = [bcol,'F_MISS',args.pheno_col]+sex_missing_cols
    pfrm.to_csv(fname,columns=cols)



noX = """*-subsection{Detailed sex analysis}

There were no X-chromosome SNPs and so this was not possible. Note that in Table*-ref{table:batchrep:all} it was not possible to compute the number of samples that failed the sex check.

"""

def detailedSexAnalysis(pfrm,missing_sex_columns):
    sex_fname = args.base+"_missing_and_sexcheck.csv"
    if missing_sex_columns == "---":
        g=open(sex_fname,"w")
        g.close()
        return noX
    dumpMissingSexTable(sex_fname,pfrm,missing_sex_columns)
    header = detSexHeader(pfrm,missing_sex_columns)
    tbl    = detSexGroup(pfrm,"overall",missing_sex_columns)
    if args.pheno_col != "all":
       g = pfrm.groupby(args.pheno_col)
       for grpname, gg in g:
           tbl = tbl+detSexGroup(gg,grpname,missing_sex_columns)
    return \
        det_sex_analysis%(sex_fname,args.pheno_col.replace('_',' ')) +\
        det_table%(header+(tbl,args.pheno_col))

       
def backslashify(text):
    return text.replace("*-",chr(92)).replace("##",chr(36))

def getBatchAnalysis():
# Show the missingness by batch
   debug("In getBatchAnalysis") #DBG
   if "${params.batch}" in no_response:
       print("Making fake batch file",ifrm.columns," index is ",ifrm.index.names)  # DBG
       args.batch_col = 'all'
       bfrm = DataFrame([1]*len(ifrm),index=ifrm.index,columns=['all'])
       res_text = "Table *-ref{table:batchrep:all} on page *-pageref{table:batchrep:all} shows the error rate."
       g=open("nopcs.pdf","w")
       g.close()
   else:
       bfrm = getCsvI(args.batch)
       debug("Read in batch file",bfrm.columns," index is ",bfrm.index.names) #DBG
       res_text = "Table *-ref{table:batchrep:%(bname)s} on page *-pageref{table:batchrep:%(bname)s} shows the error "\
                  " rate as shown by *-url{%(bname)s} as found in file *-url{%(fname)s}."\
                   %({'bname':args.batch_col,'fname':args.batch})
   
   problems = ifrm.index.difference(bfrm.index).to_series().values
   problems = list(filter(lambda fn: "_replicate_" not in fn[0] and (not fn[1].startswith("_MSK_")), problems))
   if len(problems)>0:
        print("Problem there IDs in the genotype data that are not in the batch data")
        print("You need to fix one way or the other")
        print("The IDs are")
        for p in problems:
            print(p)
        sys.exit(-1)
   result = miss_vals(ifrm,bfrm,args.batch_col,args.sexcheck_report)
   return bfrm, res_text+showResult(args.batch_col,*result)


def getPhenoAnalysis():
    pfrm = got_frame = False
    res_text = ""
    if  "${params.phenotype}" in no_response:
        if "${params.batch}" not in no_response:
            debug("Making fake phenotype file") #DBG
            args.pheno_col = 'all'
            pfrm = DataFrame(["1"]*len(ifrm),index=ifrm.index,columns=['all'])
            got_frame = True
            res_text = "Table *-ref{table:batchrep:all} on page *-pageref{table:batchrep:all} shows the *-textbf{overall} error rate."
    else:
        pfrm = getCsvI(args.phenotype)
        got_frame = True
        res_text = "Table *-ref{table:batchrep:%(bname)s} on page *-pageref{table:batchrep:%(bname)s} shows the error"\
                   " rate as shown by *-url{%(bname)s} as found in file *-url{%(fname)s}."\
                     %({'bname':args.pheno_col,'fname':args.phenotype})
    if got_frame:
        result = miss_vals(ifrm,pfrm,args.pheno_col,args.sexcheck_report)
        res_text = res_text + showResult(args.pheno_col,*result)
    else:
        pfrm = DataFrame(["1"]*len(ifrm),index=ifrm.index,columns=['all'])
        args.pheno_col="all"
    return pfrm, res_text



text = getHeader()+EOL+EOL
#if not(args.batch == "0" and args.phenotype == "0"):
ifrm = getCsvI(args.imiss)
debug("Read the missing file ",args.imiss)  # DBG
debug("Index is ",ifrm.index.names)  #DBG
debug("Columns are",ifrm.columns)  #DBG
no_response = [0,"0",False,"FALSE","false","False",""]

bfrm, btext = getBatchAnalysis()
pfrm, ptext = getPhenoAnalysis()

text = text + ptext + btext
pfrm = pfrm.join(bfrm,rsuffix="-b",lsuffix="",how='inner')
pfrm = pfrm.join(ifrm,how='inner')

col_names=['FID','IID']+list(map(lambda x: "PC%d"%x,range(1,21)))

if args.sx_pickle != "0":
    sxAnalysis = pd.read_pickle(args.sx_pickle)
    missing_sex_columns = list(map(str,sxAnalysis.columns.values))
    rdict = dict(zip(sxAnalysis.columns.values,missing_sex_columns))
    sxAnalysis.rename(columns=rdict,inplace=True)
    pfrm = pfrm.join(sxAnalysis,how='inner')
else:
    missing_sex_columns = "---"

text = text + detailedSexAnalysis(pfrm,missing_sex_columns)



eigs=getCsvI(args.eigenvec,names=col_names)

m = re.search("(.*).eigenvec",args.eigenvec)
base = m.group(1)
if null_file not in args.batch:
    figs = plotPCs(base,eigs,pfrm,args.pheno_col,bfrm,args.batch_col)
    figs_text = showFigs(figs)
else:
    figs_text=""

related_text = getRelatedPairs(pfrm,args.pheno_col,args.genome)
text = backslashify(text+related_text+figs_text)
g=open("%s-batch.tex"%args.base,"w")
g.write(text)
g.close()
   

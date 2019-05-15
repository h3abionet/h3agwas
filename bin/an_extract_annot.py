#!/usr/bin/env python3
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import sys
import pandas as pd
import argparse
import numpy as np

headerpdf="""
*-documentclass[11pt]{article}

*-usepackage[paper=a4paper,left=2cm,right=2cm,top=2cm,bottom=2cm]{geometry}
*-usepackage{graphicx}
*-usepackage{subfig}
*-usepackage{listings}
*-usepackage{longtable}
*-usepackage{array}
*-usepackage{booktabs}
*-usepackage{float}
*-usepackage{dcolumn}
*-floatstyle{ruled}
*-restylefloat{figure}
*-restylefloat{table}
*-newcommand{*-lefttblcol}{*-raggedright*-hspace{0pt}}
*-newcommand{*-righttblcol}{*-raggedleft*-hspace{0pt}}
*-newcommand{*-centretblcol}{*-centering*-hspace{0pt}}

*-newcolumntype{P}[1]{>{*-lefttblcol}p{#1}}
*-newcolumntype{Q}[1]{>{*-righttblcol}p{#1}}
*-newcolumntype{R}[1]{>{*-centretblcol}p{#1}}
*-lstset{
basicstyle=*-small*-ttfamily,
columns=flexible,
breaklines=true
}
*-usepackage{url}
*-begin{document}
"""
def GetAnnotation(File,infopos, out):
   read_info=open(File)
   out=open(out, 'w')
   head=read_info.readline()
   out.write(head)
   balisefind=False
   for Line in read_info :
     Spl=Line.split()
     if infopos[1]==Spl[0] and infopos[2]==Spl[1] :
        out.write(Line)
        balisefind=True
        break
   if balisefind==False :
      return ("",False,"")
   return (Line.replace('\n','').split(),balisefind, head.replace('\n','').split())

def GetFileForAnnot(list_file_annot):
  read_listanno=open(list_file_annot)
  File=None
  for x in read_listanno :
    tmpx=x.replace('\n','').split()
    if tmpx[0].upper()=="ALL" or tmpx[0]==infopos[1] :
      File=tmpx[1]
      break
  return File

def GetInfoHead(File) :
   DicFreq={}
   DicInfo={}
   listrs={}
   Lire=open(File)
   for Line in Lire :
     SplLin=Line.replace('\n','').split("\t")
     Type=SplLin[1].upper()
     Head=SplLin[0]
     if Type[:2]=="NA":
       continue 
     elif Type[:2]=='RS' :
       listrs[Head]=SplLin[2]
     elif Type[:2]=="IN" :
        DicInfo[Head]=SplLin[2]
     elif Type[:2]=="FR":
        SplInf=SplLin[2].split(" ")
        DB=SplInf[0]
        Type=SplInf[1].upper()
        if DB not in DicFreq :
           DicFreq[DB]={}
        if Type=="FILTER":
           DicFreq[DB]['FILTER']=SplInf[0] 
        else :
           if "FREQ" not in DicFreq[DB]:
              DicFreq[DB]["FREQ"]={}
           Pop=SplInf[2]
           if Pop not in DicFreq[DB]['FREQ'] :
              DicFreq[DB]['FREQ'][Pop]=[None,None,None]
           if SplInf[1]=="F" :
              DicFreq[DB]['FREQ'][Pop][0]=SplLin[0]
           if SplInf[1]=="N" :
              DicFreq[DB]['FREQ'][Pop][1]=SplLin[0]
           if SplInf[1]=="H" :
              DicFreq[DB]['FREQ'][Pop][2]=SplLin[0]
   return (listrs, DicFreq, DicInfo)   
           


EOL=chr(10)

def parseArguments():
    parser = argparse.ArgumentParser(description='extract annotation for specific position')
    parser.add_argument('--list_file_annot',type=str,required=True, help="file contains chro and list of files with annotation")
    parser.add_argument('--info_file_annot',type=str,required=False, help="file contains chro and list of files with annotation", default=None)
    parser.add_argument('--info_pos',type=str,required=True,help="file contain one position : rs chro pos")
    parser.add_argument('--out', type=str,help="output")
    args = parser.parse_args()
    return args

args = parseArguments()

TAB =chr(9)

args = parseArguments()
read_infopos=open(args.info_pos)
infopos=read_infopos.readline().replace("\n","").split()
rs_initname=infopos[0]

FileAnnot=GetFileForAnnot(args.list_file_annot)
latex_tex = r"""
\clearpage
\section{annotation of %(rs_initname)s analysis}
\label{sec:%(rs_initname)s}

"""

if FileAnnot==None :
   print("doesn't find chromosome or All in "+ args.info_pos+"\n")
   sys.exit(1)
(InfoAnnotPos, balisefind, head)=GetAnnotation(FileAnnot,infopos, args.out+".annot")
print(balisefind)
out_tex=args.out+".tex"

colors=['red', 'blue']
hashd={'rs_initname':rs_initname, 'listrsknow':""}
if args.info_file_annot :
   if balisefind:
      (listrs, DicFreq, DicInfo)=GetInfoHead(args.info_file_annot)
      listrsknow=": "
      for rs in listrs :
        if rs in head and  InfoAnnotPos[head.index(rs)].upper()!="NA":
           listrsknow+=InfoAnnotPos[head.index(rs)]+" " 
      latex_tex+=r"""
\subsection{annotation}
\begin{itemize}
\item rs found :%(listrsknow)s
"""
      for Info in DicInfo :
         if Info in head and  InfoAnnotPos[head.index(Info)].upper()!="NA":  
            latex_tex+="\item "+ DicInfo[Info].replace('_','-')+": "+InfoAnnotPos[head.index(Info)].replace('_','-')+"\n"
      latex_tex+="\end{itemize}"
      NbDb=len(DicFreq)
      latex_tex+=r"""
\subsection{frequencies in other population}
      """
      for DB in DicFreq.keys() :
          ## we count number frequencies
           if 'FILTER' in DicFreq[DB] and ((DicFreq[DB]['FILTER'] in head) and InfoAnnotPos[head.index(DicFreq[DB]['FILTER'])].upper()!="NA"):
              latex_tex+="Filter found in vcf file for position was " + InfoAnnotPos[head.index(DicFreq[DB]['FILTER'])]+"\n"
           CmtFreq=0
           for pop in DicFreq[DB]['FREQ'].keys() :
             headfreq=DicFreq[DB]['FREQ'][pop][0]
             if  headfreq!=None and ((headfreq in head) and InfoAnnotPos[head.index(headfreq)].upper()!="NA"):
                CmtFreq+=1 
           if CmtFreq==0 :
             continue
           nbcol=4
           numline=1
           #PosFig=int(CmtFreq/nbcol+0.5)*100 + nbcol*10 +1 
           nbrow=round(CmtFreq/nbcol+0.5)
           print(nbrow,nbcol,numline)
           cmtfigline=1
           plt.ioff()
           for pop in DicFreq[DB]['FREQ'].keys() :
              PosFig=numline*100+nbcol*10+cmtfigline
              headfreq=DicFreq[DB]['FREQ'][pop][0]
              if headfreq!=None and ((headfreq in head) and InfoAnnotPos[head.index(headfreq)].upper()!="NA"):
                 plt.suptitle(DB)
                 headN=DicFreq[DB]['FREQ'][pop][1]
                 title=DB+": "+pop
                 plt.subplot(nbrow,nbcol,cmtfigline)
                 if headN and InfoAnnotPos[head.index(headN)].upper()!="NA":
                    title+="(N: "+str(int(float(InfoAnnotPos[head.index(headN)])))+")"
                 plt.title(title, fontsize=10)
                 pi=float(InfoAnnotPos[head.index(headfreq)])
                 sizes=[pi, 1-pi]
                 plt.pie(sizes, autopct='%1.1f%%', colors=colors,shadow=True,startangle=90)
                 cmtfigline+=1
           if cmtfigline>1:
              plt.savefig(DB+"-"+rs_initname+'.pdf') 
              plt.close()
              figlatex=r"""
\subsubsection{%(DB)s}
\begin{figure}[ht]
\begin{center}
\includegraphics[width=15cm]{%(DB)s-%(rs_initname)s}
\caption{%(DB)s frequencies : frequencies for %(rs_initname)s}
\label{fig:%(DB)s-%(rs_initname)s}
\end{center}
\end{figure}
           """ 
           hashd={'DB':DB,'rs_initname':rs_initname}
           latex_tex+=figlatex%hashd
      hashd={'rs_initname':rs_initname, 'listrsknow':listrsknow}

latex_tex=latex_tex%hashd
latex_tex=headerpdf+latex_tex+"\n*-end{document}"
latex_tex=latex_tex.replace("*-",chr(92)).replace("##",chr(36)).replace("@.@",chr(10))

g = open(out_tex,"w")
g.write(latex_tex)
g.close()



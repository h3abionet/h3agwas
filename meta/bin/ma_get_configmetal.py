#!/usr/bin/env python3

import sys
import os
import argparse



def GetSepMetal(Sep):
   #WHITESPACE|COMMA|BOTH|TAB
   ListOfSep=["WHITESPACE","COMMA","TAB"]
   if len(Sep)>2 :
      Sep=Sep.upper()[:3]
      if Sep=='COM' :
         Sep='COMMA'
      elif Sep=='TAB' :
         Sep='TAB'
      elif Sep=='WHI' :
         Sep='WHITESPACE'
      else :
        return None
   else :
     if Sep==',' :
        Sep='COMMA'
     elif Sep==' ' :
        Sep='WHITESPACE'
     elif Sep==' ' :
        Sep='TAB'
   if Sep not in ListOfSep :
      return None
   return Sep

def checknull(x):
   if x=="" :
      return "NA"
   return x



def parseArguments():
    parser = argparse.ArgumentParser(description='transform file and header')
    parser.add_argument('--filelist',type=str,required=True, help="input file association")
    parser.add_argument('--output_configmetal',type=str,required=True, help="input file association")
    parser.add_argument('--out_file_metal',type=str,default=".",help="output dir if change format")
    parser.add_argument("--genomic_control", help="genomic_control ON/OFF or values see metal manual", type=str, default='F')
    parser.add_argument("--weigthedz", help="if you want stderr and not based on sample size ", type=str, default='F')
    parser.add_argument("--overlap", help="if you want Inverse Variance Weighted Meta-analysis, y or n", type=str, default='F')
    args = parser.parse_args()
    return args

#METAL_config=["MARKERLABEL", "", "","", "","EFFECTLABEL","STDERRLABEL","PVAL", "WEIGHTLABEL","FREQLABEL","STRANDLABEL","SEPARATOR"]
#head_config_i=["rsID","Chro","Pos","A1","A2","Beta","Se","Pval","N","freqA1","direction","Imputed","Sep","File"]
METAL_config=["MARKERLABEL","EFFECTLABEL","STDERRLABEL","PVAL", "WEIGHTLABEL","FREQLABEL","STRANDLABEL"]
head_config_i=["rsID","Beta","Se","Pval","N","freqA1","direction","Imputed"]


args = parseArguments()

read_listfile=open(args.filelist)
write_metal=open(args.output_configmetal, 'w')
param_MetalI=[]


if args.overlap[0] == "T":
  param_MetalI.append("##overlap\nOVERLAP ON")
  param_MetalI.append("SCHEME SAMPLESIZE")
else :
 if args.weigthedz[0]=='T' :
   param_MetalI.append("SCHEME SAMPLESIZE")
 else :
   param_MetalI.append("SCHEME STDERR")

if args.genomic_control[0]=="T" :
      param_MetalI.append("GENOMICCONTROL ON")


param_Metal=[]
nbfreq=0
for line in read_listfile :
    line=line.replace('\n', '')
    read_file=open(line)
    head_spl=[x.upper() for x in read_file.readline().replace('\n','').split()]
    read_file=read_file.close()
    for Cmt in range(len(METAL_config)) :
       if head_config_i[Cmt].upper() in head_spl :
          param_Metal.append(METAL_config[Cmt]+" "+head_config_i[Cmt]+"\n" )
    param_Metal.append("SEP TAB")
    if "DIRECTION" in head_spl:
       param_Metal.append("USESTRAND ON")
    else :
         param_Metal.append("USESTRAND OFF")
    if "A1" in head_spl and "A2" in head_spl:
         param_Metal.append("ALLELELABELS A1 A2")
    if "FREQA1" in head_spl :
        nbfreq+=1
        param_Metal.append("FREQLABEL FREQA1")
    param_Metal.append("")
    param_Metal.append("PROCESS "+line)
    param_Metal.append("")
    param_Metal.append("")


if nbfreq > 0 :
  param_MetalI.append("AVERAGEFREQ ON")
  param_MetalI.append("MINMAXFREQ ON")

param_Metal=param_MetalI+param_Metal
param_Metal+=["",""]



if args.weigthedz!='o':
     param_Metal.append("OUTFILE "+args.out_file_metal+"  .stat")
     param_Metal.append("ANALYZE\n")

param_Metal.append("QUIT\n")

write_metal=open(args.output_configmetal, 'w')
write_metal.write("\n".join(param_Metal))
write_metal.close()





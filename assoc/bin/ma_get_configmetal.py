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
    parser.add_argument("--inv_var_weigth", help="if you want Inverse Variance Weighted Meta-analysis, y or n", type=str, default='F')
    args = parser.parse_args()
    return args

#METAL_config=["MARKERLABEL", "", "","", "","EFFECTLABEL","STDERRLABEL","PVAL", "WEIGHTLABEL","FREQLABEL","STRANDLABEL","SEPARATOR"]
#head_config_i=["rsID","Chro","Pos","A1","A2","Beta","Se","Pval","N","freqA1","direction","Imputed","Sep","File"]
METAL_config=["MARKERLABEL","EFFECTLABEL","STDERRLABEL","PVAL", "WEIGHTLABEL","FREQLABEL","STRANDLABEL"]
head_config_i=["rsID","Beta","Se","Pval","N","freqA1","direction","Imputed"]


args = parseArguments()

read_listfile=open(args.filelist)
write_metal=open(args.output_configmetal, 'w')
param_Metal=[]
for line in read_listfile :
    line=line.replace('\n', '')
    read_file=open(line)
    head_spl=[x.upper() for x in read_file.readline().replace('\n','').split()]
    read_file=read_file.close()
    for Cmt in range(len(METAL_config)) :
       if head_config_i[Cmt].upper() in head_spl :
          param_Metal.append(METAL_config[Cmt]+" "+head_config_i[Cmt]+"\n" )
    param_Metal.append("SEP TAB")
    if args.inv_var_weigth[0]=='T' :
       param_Metal.append("SCHEME STDERR")
    if "DIRECTION" in head_spl:
       param_Metal.append("USESTRAND ON")
    else :
         param_Metal.append("USESTRAND OFF")
    if "A1" in head_spl and "A2" in head_spl:
         param_Metal.append("ALLELELABELS A1 A2")
    param_Metal.append("")
    param_Metal.append("PROCESS "+line)
    param_Metal.append("")
    param_Metal.append("")


param_Metal+=["",""]
if args.genomic_control[0]=="T" :
      param_Metal.append("GENOMICCONTROL " + args.genomic_control)
if args.inv_var_weigth!='o':
     param_Metal.append("OUTFILE "+args.out_file_metal+"  .stat")
     param_Metal.append("ANALYZE\n")
param_Metal.append("QUIT\n")

write_metal=open(args.output_configmetal, 'w')
write_metal.write("\n".join(param_Metal))
write_metal.close()





#!/usr/bin/env python3
''' 
manage output error for cojo in gcta
'''
from  subprocess import CalledProcessError
import subprocess
import os
import argparse
import numpy as np
import pandas as pd
import sys

EOL = chr(10)
#SNP A1 A2 freq b se p N 

def parseArguments():
    parser = argparse.ArgumentParser(description='format file for gcta, append N and frequencie if not present using bed file')
    parser.add_argument('--file_err',type=str,required=True, help="file error")
    args = parser.parse_args()
    return args


##Error: residual variance is out of boundary, the model is over-fitting. Please specify a more stringent p cutoff value.

args = parseArguments()

read_err=open(args.file_err)
TypeExit=0
for line in read_err :
    print(line)
    if "Error: residual variance" in line :
        TypeExit=14
     
sys.exit(TypeExit)


#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
@Author: Li Fajin
@Date: 2019-08-13 15:07:14
@LastEditors: Li Fajin
@LastEditTime: 2019-08-30 14:54:19
@Description: The scripts is used for merging ribosome density file generated by MetageneAnalysis.py
usage: python MergeSampleDensity.py inputFile outputFile
1) inputFile: input files separated by comma. e.g. test1_dataframe.txt,test2_dataframe.txt
2) outputFile: an output file after merging.
'''
import os
import sys
import numpy as np
import pandas as pd

inputFile=sys.argv[1]
outputFile=sys.argv[2]

def MergeSampleData(inputFile,outputFile):
    inputFiles=inputFile.strip().split(',')
    data=[pd.read_csv(File,sep="\t",header=0) for File in inputFiles]
    Final_data=pd.concat(data,axis=0)
    Final_data.to_csv(outputFile,sep="\t",index=0)

def main():
    print("Start merging...",file=sys.stderr)
    MergeSampleData(inputFile,outputFile)
    print("Finish the files merging!",file=sys.stderr)
if __name__=="__main__":
    main()
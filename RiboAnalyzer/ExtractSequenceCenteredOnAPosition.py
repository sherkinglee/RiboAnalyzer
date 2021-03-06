#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
@Author: Li Fajin
@Date: 2019-08-21 16:36:15
@LastEditors: Li Fajin
@LastEditTime: 2019-08-30 17:01:22
@Description: This script is used for extract sequence based on a specific position, for example, extract around 5 base centered on codon 50.
For example: 12345N12345
Usage: python ExtractSequenceCenteredOnAPosition.py -i cds_sequence.fa -o output_prefix --center 50 (nt) --stretch 20 (nt)
'''

import numpy as np
import sys
import os
import re
from itertools import groupby
from optparse import OptionParser
import pandas as pd
from .__init__ import __version__

def create_parser():
	usage="python %prog [options]" +"\n" + __doc__+"\n"
	parser=OptionParser(usage=usage,version=__version__)
	parser.add_option("-i",'--fasta-file',action="store",type='string',dest="cdsSequence",
		help="CDS sequences with fasta format. All sequences are from the longest transcript which may be generated by GetProteinCodingSequence.py")
	parser.add_option("--center",action="store",type="int",dest="center_position",default=150,
		help="The center position that we want to fix up. default=%default")
	parser.add_option("--stretch",action="store",type='int',dest="stretch_length",default=3,
		help="The stretch length based on the center position.default=%default.")
	parser.add_option("-o","--output_prefix",action="store",type="string",dest="output_prefix",
		help="Prefix of output files.[required]")
	return parser

def fastaIter(transcriptFile):
	fastaDict={}
	f=open(transcriptFile,'r')
	faiter=(x[1] for x in groupby(f,lambda line: line.strip()[0]==">")) ## groupby returns a tuple (key, group)
	for header in faiter:
		geneName=header.__next__().strip(">").split(" ")[0]
		seq=''.join(s.strip() for s in faiter.__next__())
		fastaDict[geneName]=seq
	return fastaDict


def extract_motif_sequence(transcriptFile,in_selectTrans,center_position,stretch_length,output_prefix):
	tranSeqDict=fastaIter(transcriptFile)
	fout=open(output_prefix+"_centered_sequence.fa",'w')
	for trans in in_selectTrans:
		cds_seq=tranSeqDict[trans][:-3] ## excluding the stop codon.
		if len(cds_seq)%3 != 0:
			continue
		if (center_position+stretch_length) > len(cds_seq):
			raise IOError("--center + --stretch <= length of a cds sequence. Please reset your parameters!")
		tmp_seq=cds_seq[int(center_position-stretch_length-1):int(center_position+stretch_length)]
		fout.write("%s%s\n" %(">",trans))
		fout.write("%s\n" %(str(tmp_seq)))

def main():
	parser=create_parser()
	(options,args)=parser.parse_args()
	(transcriptFile,center_position,stretch_length,output_prefix)=(options.cdsSequence,options.center_position,options.stretch_length,options.output_prefix)
	if not transcriptFile:
		raise IOError("Please enter the cds sequence containing interested transcripts which could be generated by GetProteinCodingSequence.py.")
	if stretch_length> center_position:
		raise IOError("Please the stretch length must be less than centered position. Namely --stretch <= --center")
	tranSeqDict=fastaIter(transcriptFile)
	in_selectTrans=tranSeqDict.keys()
	print("Start extract the sequence...")
	extract_motif_sequence(transcriptFile,in_selectTrans,center_position,stretch_length,output_prefix)
	print("Finish the step of extracting sequences.")


if __name__=="__main__":
	main()


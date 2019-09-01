#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
@Author: Li Fajin
@Date: 2019-08-18 09:26:59
@LastEditors: Li Fajin
@LastEditTime: 2019-09-01 16:00:39
@Description: This file is used for local cAI and global cAI calculation of each gene

Notes:
1) part of codes adapted from CAI package: Lee, B. D. (2018). Python Implementation of Codon Adaptation Index. Journal of Open Source Software, 3 (30), 905. https://doi.org/10.21105/joss.00905
2) input files contains:
	a. interested sequences with fasta format separated by comma.
	b. either reference cds sequences with fasta format.  e.g. some high expressed genes, whose cds sequences could be generated by GetProteinCodingSequence.py
3) output files contains:
	a. if given a reference sequence file, a RUSCs and weights used for CAI calculation would be generated.
	b. a local cAI values at each position would be generated.  One fasta file, one this kind of file.
	c. a meta results at each position would be generated.
	d. a global CAI values of each gene.

'''

from .FunctionDefinition import *
from collections import Counter
from collections import defaultdict
from itertools import chain
from functools import reduce
from operator import mul
import Bio.Data.CodonTable as ct
# get rid of Biopython warning
import warnings
from Bio import BiopythonWarning
warnings.simplefilter("ignore", BiopythonWarning)


def calculate_geometric_mean(values):
	length=len(values)
	return pow(reduce(mul,values),1/length)

def synonymous_codons(genetic_code_dict):
	''' code adapted from CAI'''

	# invert the genetic code dictionary to map each amino acid to its codons
	codons_for_amino_acid = defaultdict(list)
	for codon, amino_acid in genetic_code_dict.items():
		codons_for_amino_acid[amino_acid].append(codon)

	# create dictionary of synonymous codons
	# Example: {'CTT': ['CTT', 'CTG', 'CTA', 'CTC', 'TTA', 'TTG'], 'ATG': ['ATG']...}
	return {
		codon: codons_for_amino_acid[genetic_code_dict[codon]]
		for codon in genetic_code_dict.keys()
	}

def RSCU(sequences, genetic_code=1):
	''' code adapted from CAI'''
	if not isinstance(sequences, (list, tuple)):
		raise ValueError(
			"Be sure to pass a list of sequences, not a single sequence. "
			"To find the RSCU of a single sequence, pass it as a one element list."
		)

	# ensure all input sequences are divisible by three
	for sequence in sequences:
		if len(sequence) % 3 != 0:
			raise ValueError("Input sequence not divisible by three")
		if not sequence:
			raise ValueError("Input sequence cannot be empty")

	# count the number of each codon in the sequences
	codon_sequences = (
		(sequence[i : i + 3].upper() for i in range(0, len(sequence), 3))
		for sequence in sequences
	)
	codons = chain.from_iterable(
		codon_sequences
	)  # flat list of all codons (to be used for counting)
	counts = Counter(codons)
	# "if a certain codon is never used in the reference set... assign [its
	# count] a value of 0.5"
	for codon in ct.unambiguous_dna_by_id[genetic_code].forward_table:
		if counts[codon] == 0:
			counts[codon] = 0.5

	# determine the synonymous codons for the genetic code
	synonymous_codons = _synonymous_codons[genetic_code]

	# hold the result as it is being calulated
	result = {}

	# calculate RSCU values
	for codon in ct.unambiguous_dna_by_id[genetic_code].forward_table:
		result[codon] = counts[codon] / (
			(len(synonymous_codons[codon]) ** -1)
			* (sum((counts[_codon] for _codon in synonymous_codons[codon])))
		)


	return result


def relative_adaptiveness(sequences=None, RSCUs=None, genetic_code=1):
	''' code adapted from CAI'''

	# ensure user gave only and only one input
	if sum([bool(sequences), bool(RSCUs)]) != 1:
		raise TypeError("Must provide either reference sequences or RSCU dictionary")

	# calculate the RSCUs if only given sequences
	if sequences:
		RSCUs = RSCU(sequences, genetic_code=genetic_code)

	# determine the synonymous codons for the genetic code
	synonymous_codons = _synonymous_codons[genetic_code]

	# calculate the weights
	weight = {}
	for codon in RSCUs:
		weight[codon] = RSCUs[codon] / max(
			(RSCUs[_codon] for _codon in synonymous_codons[codon])
		)
	for stopCodon in ['TAG','TGA','TAA']:
		weight[stopCodon]=0
	return weight

def CAI_of_each_trans(sequence, weights=None, RSCUs=None, reference=None, genetic_code=1):
	''' codes adapted from CAI'''
	# validate user input
	if sum([bool(reference), bool(RSCUs)], bool(weights)) != 1:
		raise TypeError(
			"Must provide either reference sequences, or RSCU dictionary, or weights"
		)

	# validate sequence
	if not sequence:
		raise ValueError("Sequence cannot be empty")

	# make sure input sequence can be divided into codons. If so, split into list of codons
	if len(sequence) % 3 != 0:
		raise ValueError("Input sequence not divisible by three")
	sequence = sequence.upper()
	sequence = [sequence[i : i + 3] for i in range(0, len(sequence), 3)]

	# generate weights if not given
	if reference:
		weights = relative_adaptiveness(sequences=reference, genetic_code=genetic_code)
	elif RSCUs:
		weights = relative_adaptiveness(RSCUs=RSCUs, genetic_code=genetic_code)
	sequence_weights = []
	for codon in sequence:
		if codon not in _non_synonymous_codons[genetic_code]:
			try:
				sequence_weights.append(weights[codon])
			except KeyError:
				# ignore stop codons
				if codon in ct.unambiguous_dna_by_id[genetic_code].stop_codons:
					pass
				else:
					raise KeyError(
						"Bad weights dictionary passed: missing weight for codon "
						+ str(codon)
						+ "."
					)
	return float(calculate_geometric_mean(sequence_weights))

def global_cAI(sequenceDict,weights=None, RSCUs=None, reference=None, genetic_code=1):
	'''calculate global CAI'''
	cAI={}
	for trans in sequenceDict.keys():
		cds_seq=sequenceDict[trans][:-3] ## excluding stop codon
		cAI[trans]=CAI_of_each_trans(cds_seq, weights=weights, RSCUs=RSCUs, reference=reference, genetic_code=genetic_code)
	return cAI


def get_trans_frame_cAI(sequenceDict,upLength,downLength,weight,table=1):
	'''get the local cAI in each codon'''
	startcAI=np.zeros(int(upLength+downLength+1),dtype='float64')
	stopcAI=np.zeros(int(upLength+downLength+1),dtype='float64')
	startcAIList=[]
	stopcAIList=[]
	startPos=[]
	stopPos=[]
	cAI={}
	in_selectTrans=sequenceDict.keys()
	# print(weight,file=sys.stderr)
	for trans in in_selectTrans:
		cds_seq=sequenceDict[trans]
		if len(cds_seq) %3 !=0:
			continue
		tmpcAI=[]
		codon_seq=[cds_seq[i:i+3] for i in np.arange(0,len(cds_seq),3)]
		for codon in codon_seq:
			tmpcAI.append(weight[codon])
		(tmpStartWin,tmpStartPos)=getWindowsVector(upLength,downLength,tmpcAI,0) #start codon coor is 0 (0-based), codon level
		(tmpStopWin, tmpStopPos) =getWindowsVector(downLength,upLength,tmpcAI,(len(tmpcAI)-1))  #stop codon coor is len-1 (0-based) codon level
		startcAIList.append(tmpStartWin)
		stopcAIList.append(tmpStopWin)
		startPos.append(tmpStartPos)
		stopPos.append(tmpStopPos)
		cAI[trans]=tmpcAI
	startcAIList=np.array(startcAIList)
	startPos=np.array(startPos)
	stopcAIList=np.array(stopcAIList)
	stopPos=np.array(stopPos)
	for terms in np.arange(upLength+downLength+1):
		startcAI[terms]=np.mean(startcAIList[np.where(startPos[:,terms]==1),terms])
		stopcAI[terms] =np.mean(stopcAIList[np.where(stopPos[:,terms]==1),terms])
	return(startcAI,stopcAI,cAI)

def parse_weight_file(inputFile):
	'''The RSCU file contains two columns: codon and weight'''
	weights={}
	i=0
	with open(inputFile,'r') as f:
		for line in f:
			i+=1
			if line.strip()=="":
				continue
			if i==1:
				continue
			codon=str(line.strip().split("\t")[0])
			weight=float(line.strip().split("\t")[1])
			weights[codon]=weight
	return weights

def write_weight_file(weights,outputFile):
	''' output the RSCUs or weights '''
	with open(outputFile,"w") as f:
		f.write("%s\t%s\n" %("codons","weights"))
		for k,v in weights.items():
			f.write("%s\t%s\n" %(str(k),str(v)))


def write_trans_file_cAI_dataframe(inFastaAttr,outFile):
	data=[]
	for fasta in inFastaAttr:
		k=pd.DataFrame([fasta.fastaLegend]*len(fasta.startcAI))
		start=pd.DataFrame(fasta.startcAI)
		stop=pd.DataFrame(fasta.stopcAI)
		cAI=pd.merge(start,stop,how="left",left_index=True,right_index=True)
		cAI=pd.merge(k,cAI,how="left",left_index=True,right_index=True)
		data.append(cAI)
	temp=data[0]
	if len(data) < 1:
		raise EOFError("Empty file, there is nothing in the file.")
	if len(data) == 1:
		temp.columns=["sample","start_cAI","stop_cAI"]
		temp.to_csv(outFile,sep="\t",index=0)
	else:
		for i in np.arange(1,len(data)):
			temp=np.vstack((temp,data[i]))
		temp=pd.DataFrame(temp,columns=["sample","start_cAI","stop_cAI"])
		temp.to_csv(outFile,sep="\t",index=0)

def write_cAI_of_each_gene(inFastaAttr,outFile):
	data=[]
	data_index=[]
	for fasta in inFastaAttr:
		d=fasta.cAI
		i=fasta.fastaLegend
		data.append(d)
		data_index.append(i)
	data=pd.DataFrame(data,index=data_index)
	data=data.T
	data.to_csv(outFile,sep="\t")

def write_cAI_per_codon(inFastaAttr,outFile):
	for fasta in inFastaAttr:
		with open(outFile+"_"+fasta.fastaLegend+"_local_cAI_each_position.txt",'w') as f:
			f.write("%s\t%s\n" %("transcripts","local_cAI"))
			for trans,cAIperCodons in fasta.cAIPerCodon.items():
				f.write("%s\t" %(trans))
				for pos in range(len(cAIperCodons)):
					f.write("%s\t" %(str(cAIperCodons[pos])))
				f.write("\n")


def main():
	parser=create_parser_for_cAI()
	(options,args)=parser.parse_args()
	_synonymous_codons = {k: synonymous_codons(v.forward_table) for k, v in ct.unambiguous_dna_by_id.items()}
	_non_synonymous_codons = {k: {codon for codon in v.keys() if len(v[codon]) == 1}for k, v in _synonymous_codons.items()}
	# validate user input
	if sum([bool(options.reference), bool(options.RSCUs)], bool(options.weights)) != 1:
		raise TypeError(
			"Must provide either reference sequences, or RSCU dictionary, or weights"
		)
	if not options.transcriptFiles:
		raise IOError("Please input your interested sequences with fasta format!")
	transcriptFiles=options.transcriptFiles.strip().split(",")
	transFileLegend=options.trans_file_legend.strip().split(",")
	references=options.reference
	RSCUs=options.RSCUs
	weights=options.weights
	## handle bam file attr
	fasta_attr=[]
	for ii,jj in zip(transcriptFiles,transFileLegend):
		fasta=fasta_attrbution(ii,jj)
		fasta_attr.append(fasta)
	print("your input : "+ str(len(transcriptFiles))+" transcript files",file=sys.stderr)
	if references:
		referenceDict=fastaIter(references)
		references=list(referenceDict.values())
		RSCUs=RSCU(references,genetic_code=options.genetic_table)
		weights = relative_adaptiveness(RSCUs=RSCUs, genetic_code=options.genetic_table)
		write_weight_file(weights,options.output_prefix+"_weights_for_cAI.txt")
		write_weight_file(RSCUs,options.output_prefix+"_RSCUs_for_cAI.txt")
		for fasta in fasta_attr:
			sequenceDict=fastaIter(fasta.fastaName)
			(fasta.startcAI,fasta.stopcAI,fasta.cAIPerCodon) = get_trans_frame_cAI(sequenceDict,options.upstream_codon,options.downstream_codon,weights,table=options.genetic_table)
			fasta.cAI=global_cAI(sequenceDict,weights=weights,genetic_code=options.genetic_table)
	elif RSCUs:
		RSCUs=parse_weight_file(RSCUs)
		weights = relative_adaptiveness(RSCUs=RSCUs, genetic_code=options.genetic_table)
		write_weight_file(weights,options.output_prefix+"_weights_for_cAI.txt")
		for fasta in fasta_attr:
			sequenceDict=fastaIter(fasta.fastaName)
			(fasta.startcAI,fasta.stopcAI,fasta.cAIPerCodon) = get_trans_frame_cAI(sequenceDict,options.upstream_codon,options.downstream_codon,weights,table=options.genetic_table)
			fasta.cAI=global_cAI(sequenceDict,RSCUs=RSCUs,genetic_code=options.genetic_table)
	elif weights:
		weights=parse_weight_file(weights)
		for fasta in fasta_attr:
			sequenceDict=fastaIter(fasta.fastaName)
			(fasta.startcAI,fasta.stopcAI,fasta.cAIPerCodon) = get_trans_frame_cAI(sequenceDict,options.upstream_codon,options.downstream_codon,weights,table=options.genetic_table)
			fasta.cAI=global_cAI(sequenceDict,weights=weights,genetic_code=options.genetic_table)
	else:
		raise IOError("Please enter a correct input. Must provide either reference sequences, or RSCU dictionary, or weights")


	print("Finish the step of get_trans_frame_cAI...",file=sys.stderr)
	write_trans_file_cAI_dataframe(fasta_attr,options.output_prefix+"_local_cAI_dataframe.txt")
	write_cAI_of_each_gene(fasta_attr,options.output_prefix+"_global_cAI.txt")
	write_cAI_per_codon(fasta_attr,options.output_prefix)
	print("Finish the step of cAI calculation.",file=sys.stderr)

if  __name__=="__main__":
	main()

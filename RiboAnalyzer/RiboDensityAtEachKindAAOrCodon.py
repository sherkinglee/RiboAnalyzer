#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
@Author: Li Fajin
@Date: 2019-08-19 23:01:51
@LastEditors: Li Fajin
@LastEditTime: 2019-08-30 17:08:08
@Description: This script used for calculating ribosome density at each kind of amino acid or codon.
And the ribosome density means counts/sequence_depth.
input:
	the input are just the same as MetageneAnalysis.py but with longest.transcripts.cds.sequence.fa file input also which could be generated by GetProteinCodingSequence.py
output:
	if -u and -d not set: there is only one output file with Row x Col where Row is the numer of  codon kinds and the Col is the number of columns which usually are just like:

	AA      codon   si_ctrl_1_80S   si_ctrl_2_80S   si_ctrl_3_80S   si_3e_1_80S     si_3e_2_80S     si_3e_3_80S
	K       AAA     0.01689667937572358     0.01592644427516077     0.015870156244881917    0.013608564821264385    0.013323948113603084    0.122
	N       AAC     0.020736512248088333    0.02236358252621296     0.0213430945439855      0.02247569028447243     0.02225038863195103     0.234
	K       AAG     0.03195537918496984     0.03208825472187564     0.032316267592562266    0.02973099105625366     0.0282806341957415      0.1212
	N       AAT     0.019226280160431212    0.019885910677272484    0.020112024632316816    0.019848987211278594    0.019645264495487997    0.123

	where the first column is AA, the second column is codon, and the rest of other columns are samples.
	if -u and -d are set: there are two output files with the same format as what we talked above.
	the first one file is the density on the specific region and the other file is the density out of the specific region.
'''


from .FunctionDefinition import *
import re
from collections import defaultdict
from itertools import groupby
import Bio.Data.CodonTable as ct

def codon_density(in_bamFile,in_selectTrans,in_transLengthDict,in_startCodonCoorDict,in_stopCodonCoorDict,in_readLengths,in_readOffset,inCDS_countsFilterParma,inCDS_lengthFilterParma,transcriptFasta,left_pos,right_pos,mode):
	'''output the density of each class of codon'''
	transcript_sequence=fastaIter(transcriptFasta)
	passTransSet=set()
	filter_1=0
	filter_2=0
	filter_3=0
	if left_pos and (not right_pos):
		raise IOError("Please input the right position match with the left position!")
	if (not left_pos) and right_pos:
		raise IOError("Please input the left position match with the right position!")
	if (not left_pos) and (not right_pos):
		pysamFile=pysam.AlignmentFile(in_bamFile,"rb")
		pysamFile_trans=pysamFile.references
		in_selectTrans=set(pysamFile_trans).intersection(in_selectTrans)
		all_codon_density=defaultdict(list)
		all_counts=0
		for trans in in_selectTrans:
			leftCoor =int(in_startCodonCoorDict[trans])-1
			rightCoor=int(in_stopCodonCoorDict[trans])-3
			(trans_counts,read_counts_frameSum,total_reads,cds_reads)=get_trans_frame_counts(pysamFile, trans, in_readLengths, in_readOffset, in_transLengthDict[trans], leftCoor, rightCoor)
			all_counts+=total_reads
		for trans in in_selectTrans:
			leftCoor =int(in_startCodonCoorDict[trans])-1 #the first base of start codon 0-base
			rightCoor=int(in_stopCodonCoorDict[trans])-3 #the first base of stop codon 0-base
			CDS_seq=transcript_sequence[trans][:-3] ## excluding the stop codon
			# print(CDS_seq)
			if len(CDS_seq) % 3 != 0:
				filter_1+=1
				continue
			if len(CDS_seq) < inCDS_lengthFilterParma:
				filter_2+=1
				continue
			codon_seq=[]
			for i in np.arange(0,len(CDS_seq),3):
				codon=CDS_seq[i:(i+3)]
				codon_seq.append(codon)
			# print(codon_seq)
			(read_counts,read_counts_frameSum,trans_reads,cds_reads)=get_trans_frame_counts(pysamFile, trans, in_readLengths, in_readOffset, in_transLengthDict[trans], leftCoor, rightCoor)
			cds_reads_norm = 10**9*(cds_reads/(all_counts*len(read_counts_frameSum)*3))
			if mode == 'RPKM':
				if cds_reads_norm < inCDS_countsFilterParma:
					filter_3+=1
					continue
			if mode == 'counts':
				if cds_reads < inCDS_countsFilterParma:
					filter_3+=1
					continue
			all_read_counts_frame=np.array(read_counts_frameSum)
			codon_seq=np.array(codon_seq)
			# all_densityPerTrans=10**6*(all_read_counts_frame/all_counts)
			all_densityPerTrans=all_read_counts_frame/all_counts
			# print(trans+"\t"+"codon: "+str(len(codon_seq))+"\t"+"counts: "+str(len(all_densityPerTrans))+"\t" + "TAG: " + str(list(codon_seq).count('TAG'))+"\t" + "TGA: " + str(list(codon_seq).count('TGA'))+"\t" + "TAA: " + str(list(codon_seq).count('TAA')),file=sys.stderr)
			for k,v in zip(codon_seq,all_densityPerTrans):
				all_codon_density[k].append(v)
			passTransSet.add(trans)
		pysamFile.close()
		for k,v in all_codon_density.items():
			all_codon_density[k]=np.sum(v)
		print("The number of transcripts filtered by filter1 (len(CDS_seq)%3!=0) is : " + str(filter_1),file=sys.stderr)
		print("The number of transcripts filtered by filter2 (len(CDS_seq)<inCDS_lengthFilterParma) is :" + str(filter_2),file=sys.stderr)
		print("The number of transcripts filtered by filter3 (cds_reads < inCDS_countsFilterParma) is : "+str(filter_3),file=sys.stderr)
		print("The number of transcripts used for following analysis is: " + str(len(passTransSet)),file=sys.stderr)
		return all_codon_density
	if left_pos and right_pos:
		if right_pos <= left_pos:
			raise IOError("The right position must be larger than the left position!")
		pysamFile=pysam.AlignmentFile(in_bamFile,"rb")
		pysamFile_trans=pysamFile.references
		in_selectTrans=set(pysamFile_trans).intersection(in_selectTrans)
		specific_range_codon_density = defaultdict(list)
		the_rest_range_codon_density=defaultdict(list)
		all_counts=0
		for trans in in_selectTrans:
			leftCoor =int(in_startCodonCoorDict[trans])-1
			rightCoor=int(in_stopCodonCoorDict[trans])-3
			(trans_counts,read_counts_frameSum,total_reads,cds_reads)=get_trans_frame_counts(pysamFile, trans, in_readLengths, in_readOffset, in_transLengthDict[trans], leftCoor, rightCoor)
			all_counts+=total_reads
		for trans in in_selectTrans:
			if (int(right_pos)-int(left_pos)) > in_transLengthDict[trans]:
				raise IOError("The length of the range you chose is larger than some transcript! please reset your parameters!")
			if right_pos > in_transLengthDict[trans]:
				raise IOError("The right position is larger than the length of some transcrip! please reset your parameters!")
			leftCoor =int(in_startCodonCoorDict[trans])-1 #the first base of start codon 0-base
			rightCoor=int(in_stopCodonCoorDict[trans])-3 #the first base of stop codon 0-base
			CDS_seq=transcript_sequence[trans][:-3]
			if len(CDS_seq) % 3 != 0:
				filter_1+=1
				continue
			if len(CDS_seq) < inCDS_lengthFilterParma:
				filter_2+=1
				continue
			codon_seq=[]
			for i in np.arange(0,len(CDS_seq),3):
				codon=CDS_seq[i:(i+3)]
				codon_seq.append(codon)
			(read_counts,read_counts_frameSum,trans_reads,cds_reads)=get_trans_frame_counts(pysamFile, trans, in_readLengths, in_readOffset, in_transLengthDict[trans], leftCoor, rightCoor)
			cds_reads_norm = 10**9*(cds_reads/(all_counts*len(read_counts_frameSum)*3))
			if mode == 'RPKM':
				if cds_reads_norm < inCDS_countsFilterParma:
					filter_3+=1
					continue
			if mode == 'counts':
				if cds_reads < inCDS_countsFilterParma:
					filter_3+=1
					continue
			specific_range_read_counts_frame=np.array(read_counts_frameSum[int(left_pos-1):right_pos])
			the_rest_range_read_counts_frame=np.array(read_counts_frameSum[right_pos:])
			codon_seq=np.array(codon_seq)
			specific_range_codon_seq=codon_seq[int(left_pos-1):right_pos]
			the_rest_range_codon_seq=codon_seq[right_pos:]
			specific_range_densityPerTrans=specific_range_read_counts_frame/all_counts
			the_rest_range_densityPerTrans=the_rest_range_read_counts_frame/all_counts
			# specific_range_densityPerTrans=10**6*(specific_range_read_counts_frame/all_counts)
			# the_rest_range_densityPerTrans=10**6*(the_rest_range_read_counts_frame/all_counts)
			for k,v in zip(specific_range_codon_seq,specific_range_densityPerTrans):
				specific_range_codon_density[k].append(v)
			for k,v in zip(the_rest_range_codon_seq,the_rest_range_densityPerTrans):
				the_rest_range_codon_density[k].append(v)
			passTransSet.add(trans)
		pysamFile.close()
		for k,v in specific_range_codon_density.items():
			specific_range_codon_density[k]=np.sum(v)
		for k,v in the_rest_range_codon_density.items():
			the_rest_range_codon_density[k]=np.sum(v)
		print("The number of transcripts filtered by filter1 (len(CDS_seq)%3!=0) is : " + str(filter_1),file=sys.stderr)
		print("The number of transcripts filtered by filter2 (len(CDS_seq)<inCDS_lengthFilterParma) is :" + str(filter_2),file=sys.stderr)
		print("The number of transcripts filtered by filter3 (cds_reads < inCDS_countsFilterParma) is : "+str(filter_3),file=sys.stderr)
		print("The number of transcripts used for following analysis is: " + str(len(passTransSet)),file=sys.stderr)
		return specific_range_codon_density,the_rest_range_codon_density

def codon2AA(table=1):
	''' table=1 is a standard codon table'''
	codon2AA_dict={k: v.forward_table for k, v in ct.unambiguous_dna_by_id.items()}
	return codon2AA_dict[table]
def get_stop_codons(table=1):
	stop_codons_dict={k: v.stop_codons for k, v in ct.unambiguous_dna_by_id.items()}
	return stop_codons_dict[table]
def shapeData(data,table=1):
	''' use standard codon table'''
	Codon_table=codon2AA(table=table)
	stopCodons=get_stop_codons(table=table)
	for stop_codon in stopCodons:
		Codon_table[stop_codon]="*"
	data=data.reset_index()
	data=data.rename(columns={'index':'codon'})
	codon_list=data.iloc[:,0].values
	AA_list=[]
	for codon in codon_list:
		AA_list.append(Codon_table[codon])
	AA_list=pd.DataFrame(AA_list,columns=['AA'])
	data=pd.concat((AA_list,data),axis=1)
	return data


def write_specific_range_codon_density(inBamAttr,outFile,table=1):
	specific_data=[]
	rest_data=[]
	data_index=[]
	for bms in inBamAttr:
		specific_codon_density=bms.specific_range_codon_density
		rest_codon_density=bms.the_rest_range_codon_density
		sample=bms.bamLegend
		specific_data.append(specific_codon_density)
		rest_data.append(rest_codon_density)
		data_index.append(sample)
	specific_data=pd.DataFrame(specific_data,index=data_index)
	rest_data=pd.DataFrame(rest_data,index=data_index)
	specific_data=specific_data.T
	rest_data=rest_data.T
	specific_data=shapeData(specific_data,table=table)
	rest_data=shapeData(rest_data,table=table)
	specific_data.to_csv(outFile+"_specific_range_codon_density.txt",sep='\t',index=0)
	rest_data.to_csv(outFile+"_rest_range_codon_density.txt",sep="\t",index=0)

def write_all_codon_density(inBamAttr,outFile,table=1):
	all_data=[]
	data_index=[]
	for bms in inBamAttr:
		all_codon_density=bms.all_codon_density
		sample=bms.bamLegend
		all_data.append(all_codon_density)
		data_index.append(sample)
	all_data=pd.DataFrame(all_data,index=data_index)
	all_data=all_data.T
	all_data=shapeData(all_data,table=table)
	all_data.to_csv(outFile+"_all_codon_density.txt",sep="\t",index=0)

def parse_args_for_codon_density_calculation():
	parsed=creat_parser_for_calculation_of_codon_density()
	(options,args)=parsed.parse_args()
	(output_prefix, min_cds_codon,min_cds_counts,left_position, right_position, transcript_fasta,mode,table) = (options.output_prefix,options.min_cds_codon,
	options.min_cds_counts,options.upstream_position,options.downstream_position,options.transcript_fasta,options.mode,options.geneticCode)
	if options.bamListFile and (options.bam_files or options.read_length or options.read_offset or options.bam_file_legend):
		raise IOError("'-f' parameter and '-i -r -s -t' are mutually exclusive.")
	if options.bamListFile:
		bamFiles,readLengths,Offsets,bamLegends=parse_bamListFile(options.bamListFile)
	elif options.bam_files:
		bamFiles,readLengths,Offsets,bamLegends=options.bam_files.split(","),options.read_length.split("_"),options.read_offset.split("_"),options.bam_file_legend.split(",")
	else:
		raise IOError("Please check you input files!")
	print("your input : "+ str(len(bamFiles))+" bam files",file=sys.stderr)
	bam_attr=[]
	for ii,jj,mm,nn in zip(bamFiles,readLengths,Offsets,bamLegends):
		bam=bam_file_attr(ii,jj,mm,nn)
		bam_attr.append(bam)
	selectTrans,transLengthDict,startCodonCoorDict,stopCodonCoorDict,transID2geneID,transID2geneName,cdsLengthDict=reload_transcripts_information(options.coorFile)
	geneID2transID={v:k for k,v in transID2geneID.items()}
	geneName2transID={v:k for k,v in transID2geneName.items()}
	if options.in_selectTrans:
		select_trans=pd.read_csv(options.in_selectTrans,sep="\t")
		select_trans=set(select_trans.iloc[:,0].values)
		if options.id_type == 'transcript_id':
			select_trans=select_trans.intersection(selectTrans)
			print("There are " + str(len(select_trans)) + " transcripts from "+options.in_selectTrans+" used for following analysis.",file=sys.stderr)
		elif options.id_type == 'gene_id':
			tmp=[geneID2transID[gene_id] for gene_id in select_trans]
			select_trans=set(tmp)
			select_trans=select_trans.intersection(selectTrans)
			print("There are " + str(len(select_trans))+" gene id could be transformed into transcript id and used for following analysis.",file=sys.stderr)
		elif options.id_type == 'gene_name' or options.id_type=='gene_symbol':
			tmp=[geneName2transID[gene_name] for gene_name in select_trans]
			select_trans=set(tmp)
			select_trans=select_trans.intersection(selectTrans)
			print("There are " + str(len(select_trans))+" gene symbol could be transformed into transcript id and used for following analysis.",file=sys.stderr)
		else:
			raise IOError("Please input a approproate id_type parameters.[transcript_id/gene_id/gene_name/]")
	else:
		select_trans=selectTrans

	if left_position and right_position:
		for bamfs in bam_attr:
			(bamfs.specific_range_codon_density,bamfs.the_rest_range_codon_density)=codon_density(bamfs.bamName,select_trans,transLengthDict,startCodonCoorDict,
			stopCodonCoorDict,bamfs.bamLen,bamfs.bamOffset,min_cds_counts,min_cds_codon,transcript_fasta,left_position,right_position,mode)
		write_specific_range_codon_density(bam_attr,output_prefix,table=table)
		print("Finish the step of codon density calculation",file=sys.stderr)
	if (not left_position) and (not right_position):
		for bamfs in bam_attr:
			bamfs.all_codon_density=codon_density(bamfs.bamName,select_trans,transLengthDict,startCodonCoorDict,
			stopCodonCoorDict,bamfs.bamLen,bamfs.bamOffset,min_cds_counts,min_cds_codon,transcript_fasta,left_position,right_position,mode)
		write_all_codon_density(bam_attr,output_prefix,table=table)
		print("Finish the step of codon density calculation",file=sys.stderr)


def main():
	"""main program"""
	parse_args_for_codon_density_calculation()

if __name__ == "__main__":
		main()

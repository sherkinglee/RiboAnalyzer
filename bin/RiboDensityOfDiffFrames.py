#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
@Author: Li Fajin
@Date: 2019-08-16 08:51:14
@LastEditors: Li Fajin
@LastEditTime: 2019-08-16 14:33:55
@Description: This script is used for calculating ribosome density for each different reading frame.
'''


from FunctionDefinition import *
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from optparse import OptionParser

def create_parser_for_frame_density():
	'''argument parser.'''
	usage="usage: python %prog [options]" + '\n' + __doc__ + "\n"
	parser=OptionParser(usage=usage,version="%prog")
	parser.add_option("-f","--bamListFile",action="store",type="string",default=None,dest="bamListFile",
			help="Bam file list, containing 4 columns.Namely bamFiles,readLength, offSet, bamLegend. '-f' and '-i, -r, -s, -t' parameters are mutually exclusive.default=%default.")
	parser.add_option("-i","--input", action="store",type="string",default=None,dest="bam_files",
			help="Input file(s) in bam format. All files should be split by comma e.g. 1.bam,2.bam,3.bam[required]. '-i' and '-f' are mutually exclusive. default=%default")
	parser.add_option("-c","--coordinateFile",action="store",type="string",dest="coorFile",
			help="The file should contain the coordinate of start and stop codon. Generated by OutputTranscriptInfo.py.[required]")
	parser.add_option("-o","--otput_prefix",action="store",type="string",dest="output_prefix",
			help="Prefix of output files.[required]")
	parser.add_option("-r","--specific_reads_length",action="store",type="string",dest="read_length",
			help="Specific the lenght to do analysis, comma split. e.g. '28,29,30'.If use all length set 'All'. Bam files diff length select split by '_' e.g. '28,29,30_ALL_27,28' [required]. '-r' and '-f' are mutually exclusive.")
	parser.add_option("-s","--offset",action="store",type="string",dest="read_offset",
			help="Specific the offset corresponding to read length, comma split. e.g. '12,13,13'. No offset set 0. Bam files diff offset select split by '_' e.g. '12,13,13_0_12,12' [required]. '-s' and '-f' are mutually exclusive.")
	parser.add_option("-t","--bam_file_legend",action="store",type="string",dest="bam_file_legend",
			help="The legend of each bam files, comma split. e.g. 'condition1,condition2,condition3' [required]. '-t' and '-f' are mutually exclusive.")
	parser.add_option('-S','--select_trans_list',action="store",type='string',dest='in_selectTrans',
			help="Selected transcript list used for metagene analysis.This files requires the first column must be the transcript ID  with a column name.")
	parser.add_option('--id-type',action="store",type="string",dest="id_type",default="transcript_id",
			help="define the id type users input. the default is transcript id, if not, will be transformed into transcript id. default=%default")
	parser.add_option("--plot",action="store",type="string",dest="plot",default='yes',
			help="Output a  plot or not.default=%default")
	return parser



def get_ribo_density_of_different_frame(ribo_fileobj, transcript_name, read_lengths, read_offsets, transLength, startCoor, stopCoor):
	"""For each mapped read of the given transcript in the BAM file,get the P-site and codon unit reads density
	ribo_fileobj -- file object - BAM file opened using pysam AlignmentFile
	transcript_name -- Name of transcript to get counts for
	read_length -- If provided, get counts only for reads of this length.
	read_offsets -- the offset length corresponding to 5' mapped position.
	transLength -- the length of the transcript.
	startCoor -- the coordinate of the first base of start codon 0-based.
	stopCoor -- the coordinate of the first base of stop codon 0-based.
	"""
	read_counts = np.zeros(transLength,dtype="int64")
	total_reads = 0
	if read_lengths == "ALL" : ## RNA
		for record in ribo_fileobj.fetch(transcript_name):
			if record.flag == 16 or record.flag == 272:
				continue
			total_reads += 1
			position = record.pos
			read_counts[position]+=1
	else:
		read_lengths=lengths_offsets_split(read_lengths)
		read_offsets=lengths_offsets_split(read_offsets)
		for record in ribo_fileobj.fetch(transcript_name):
			if record.flag == 16 or record.flag == 272:
				continue
			for R_length, R_offset in zip(read_lengths,read_offsets):
				if  record.query_length == R_length :
					# if an offset is specified, increment position by that offset.
					position = record.pos + R_offset ## transform into the position of P-site
				else:
					# ignore other reads/lengths
					continue
				total_reads += 1
				try:
					read_counts[position]+=1
				except KeyError:
					print("Dont has this position after offset : transcript_name -> position"+" "+transcript_name+" -> "+position)
	#get trans counts for each 3 frames
	read_counts_frame0=read_counts[(startCoor+0):(stopCoor-2):3]
	read_counts_frame1=read_counts[(startCoor+1):(stopCoor-1):3]
	read_counts_frame2=read_counts[(startCoor+2):(stopCoor-0):3]
	frame0=np.sum(read_counts_frame0)
	frame1=np.sum(read_counts_frame1)
	frame2=np.sum(read_counts_frame2)
	frameSum=np.sum(frame0+frame1+frame2)
	return frame0,frame1,frame2,frameSum

def Output_frame_density(in_bamFile,in_selectTrans,in_transLengthDict,in_startCodonCoorDict,in_stopCodonCoorDict,in_readLengths,in_readOffset,output_prefix):
	'''
	output density on different frames.
	echo row means each transcripts and four columns represents frame0,frame1,frame2 and framSum.
	'''
	pysamFile=pysam.AlignmentFile(in_bamFile,"rb")
	frame_density=defaultdict(int)
	with open(output_prefix+"_reading_frames.txt",'w') as f:
		f.write("%s\t%s\t%s\t%s\t%s\n" %("transcript_id","frame0","frame1","frame2","frameSum"))
		for trans in in_selectTrans:
			leftCoor =int(in_startCodonCoorDict[trans])-1
			rightCoor=int(in_stopCodonCoorDict[trans])-3
			frame0,frame1,frame2,frameSum=get_ribo_density_of_different_frame(pysamFile, trans, in_readLengths, in_readOffset, in_transLengthDict[trans], leftCoor, rightCoor)
			f.write("%s\t%s\t%s\t%s\t%s\n" %(trans,str(frame0),str(frame1),str(frame2),str(frameSum)))
			frame_density['frame0']+=frame0
			frame_density['frame1']+=frame1
			frame_density['frame2']+=frame2

	return frame_density

def plot_bar(frame_density,output_prefix):
	''' plot density from different frame'''
	text_font={"size":20,"family":"Arial","weight":"bold"}
	plt.rc('font',weight='bold')
	fig=plt.figure(figsize=(5,4))
	ax=fig.add_subplot(111)
	colors=sns.color_palette('husl',len(frame_density))
	plt.bar(frame_density.keys(),frame_density.values(),color=colors)
	ax.set_ylabel("Raw counts",fontdict=text_font)
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.spines["bottom"].set_linewidth(2)
	ax.spines["left"].set_linewidth(2)
	ax.tick_params(which="both",width=2,length=2)
	plt.tight_layout()
	plt.savefig(output_prefix+"_reading_frames.pdf")
	plt.close()

def parse_args_for_frame_density():
	parsed=create_parser_for_frame_density()
	(options,args)=parsed.parse_args()
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
	## calculate density for each bam files
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

	for bamfs in bam_attr:
		(bamfs.frame_density) = Output_frame_density(bamfs.bamName,select_trans,transLengthDict,startCodonCoorDict,stopCodonCoorDict,bamfs.bamLen,bamfs.bamOffset,options.output_prefix+"_"+bamfs.bamLegend)
		if options.plot.upper() in ['YES','Y']:
			plot_bar(bamfs.frame_density,options.output_prefix+"_"+bamfs.bamLegend)
			print("Finish the step of plotting!",file=sys.stderr)
		elif options.plot.upper() in ['NO','N','NONE']:
			pass
		else:
			raise IOError("Please input a correct --plot parameter! [yes/no]")

	print("Finish the step of RiboDensityOfDiffFrames!",file=sys.stderr)



def main():
	"""main program"""
	## parse the options and args and output
	parse_args_for_frame_density()


if __name__ == "__main__":
		main()
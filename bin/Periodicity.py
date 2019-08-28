#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
@Author: Li Fajin
@Date: 2019-08-15 21:29:12
@LastEditors: Li Fajin
@LastEditTime: 2019-08-26 20:37:02
@Description: This script is used for checking periodicity of ribosome profiling data, but without P-site identification.
And the part code are adapted from RiboCode our lab developed before. [Xiao, et al. NAR.2018]
usage: python Periodicity -i bam -c longest.trans.info.txt -o outprefix -L 25 -R 35 --id-type transcript-id
'''


from FunctionDefinition import *
from functools import reduce
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from optparse import OptionParser


def create_parser_for_periodicity():
	'''argument parser.'''
	usage="usage: python %prog [options]" + "\n"
	parser=OptionParser(usage=usage)
	parser.add_option("-i","--input", action="store",type="string",dest="bamFile",
			help="Input file(s) in bam format.")
	parser.add_option("-c","--coordinateFile",action="store",type="string",dest="coorFile",
			help="The file should contain the coordinate of start and stop codon. Generated by OutputTranscriptInfo.py.[required]")
	parser.add_option("-o","--otput_prefix",action="store",type="string",dest="output_prefix",
			help="Prefix of output files.[required]")
	parser.add_option("-L","--left_length",action="store",type="int",default=25,dest="left_length",help="The left range of read length we will consider.")
	parser.add_option("-R","--right_length",action="store",type="int",default=35,dest="right_length",help="The right range of read length we will consider.")
	parser.add_option('-S','--select_trans_list',action="store",type='string',dest='in_selectTrans',
			help="Selected transcript list used for metagene analysis.This files requires the first column must be the transcript ID  with a column name.")
	parser.add_option('--id-type',action="store",type="string",dest="id_type",default="transcript_id",
			help="define the id type users input. the default is transcript id, if not, will be transformed into transcript id. %default=default")
	return parser

def periodicity(in_bamFile,in_selectTrans,in_transLengthDict,in_startCodonCoorDict,in_stopCodonCoorDict,left_length,right_length):
	pysamFile=pysam.AlignmentFile(in_bamFile,"rb")
	pysamFile_trans=pysamFile.references
	in_selectTrans=set(pysamFile_trans).intersection(in_selectTrans)
	start_density=defaultdict(lambda:np.zeros(101,dtype="int64"))
	stop_density=defaultdict(lambda:np.zeros(101,dtype="int64"))
	total_reads=0
	specific_counts=defaultdict(int)
	i=0
	for trans in in_selectTrans:
		i+=1
		leftCoor =int(in_startCodonCoorDict[trans])-1 #the first base of start codon 0-base
		rightCoor=int(in_stopCodonCoorDict[trans])-3 #the first base of stop codon 0-base
		# print(trans)
		for record in pysamFile.fetch(trans):
			if record.flag == 16 or record.flag == 272:
				continue
			total_reads += 1
			specific_counts[record.query_length]+=1
			if left_length<=record.query_length<=right_length:
				distance_to_start=record.reference_start-leftCoor
				dsitance_to_stop=record.reference_start-rightCoor+3
				if abs(distance_to_start)<=50:
					start_density[record.query_length][50+distance_to_start]+=1
				if abs(dsitance_to_stop) <= 50:
					stop_density[record.query_length][50+dsitance_to_stop]+=1

	specific_counts['sum']=total_reads
	print("There are " + str(i)+" used for statistics.",file=sys.stderr)
	return start_density,stop_density,total_reads,specific_counts

def flatten(xs):
	for x in xs:
		if isinstance(x,tuple):
			for xx in flatten(x):
				yield xx
		else:
			yield x

def plot_periodicity(start_density,stop_density,specific_counts,output_prefix,left_length,right_length):
	''' Part scripts adapted from RiboCode'''
	start_density=pd.DataFrame(start_density)
	stop_density=pd.DataFrame(stop_density)
	start_density=start_density.reindex(columns=sorted(start_density.columns))
	stop_density=stop_density.reindex(columns=sorted(stop_density.columns))
	start_density['sum']=start_density.apply(sum,axis=1)
	stop_density['sum']=stop_density.apply(sum,axis=1)
	start_density.loc['sum']=np.array([specific_counts[key] for key in start_density.columns])
	stop_density.loc['sum']=np.array([specific_counts[key] for key in stop_density.columns])
	start_density.to_csv(output_prefix+"_from_start_codon.txt",index=0,sep="\t")
	stop_density.to_csv(output_prefix+"_from_stop_codon.txt",index=0,sep="\t")
	with PdfPages(output_prefix + "_periodicity.pdf") as pdf:
		x = np.arange(-50,51,dtype=int)
		colors = np.tile(["b","g","r"], 34)
		for L in start_density.columns:
			xticks=[-40,-20,0,20,40]
			perct = '{:.2%}'.format(specific_counts[L] / specific_counts['sum'])
			fig,(ax1,ax2) = plt.subplots(nrows=2,ncols=1)
			y1=start_density.loc[:,L]
			y2=stop_density.loc[:,L]
			ax1.vlines(x,ymin=np.zeros(101),ymax=y1,colors=colors[:-1])
			ax1.tick_params(axis='x',which="both",top=False,direction='out')
			ax1.set_xticks(xticks)
			ax1.set_xlim((-50,50))
			ax1.set_xlabel("Distance (nt)")
			ax1.set_ylabel("Alignments")
			if not L == 'sum':
				ax1.set_title("({} nt reads,proportion:{})".format(L,perct) + "\n Distance 5'- start codons")
			else:
				ax1.set_title("Total reads, proporation:{}".format(perct) +"\n Distance 5'- start codons")
			ax2.vlines(x,ymin=np.zeros(101),ymax=y2,colors=colors[:-1])
			ax2.tick_params(axis='x',which="both",top=False,direction='out')
			ax2.set_xticks(xticks)
			ax2.set_xlim((-50,50))
			ax2.set_xlabel("Distance (nt)")
			ax2.set_ylabel("Alignments")
			ax2.set_title("Distance 5'- stop codons")

			fig.tight_layout()
			pdf.savefig(fig)
			plt.close()

	return None


def parse_args():
	parsed=create_parser_for_periodicity()
	(options,args)=parsed.parse_args()
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
	print("There are "+str(len(select_trans))+" transcripts with both start and stop codon will be used for following analysis.",file=sys.stderr)
	print("Start calculating read density...",file=sys.stderr)
	start_density,stop_density,total_reads,specific_counts=periodicity(options.bamFile,select_trans,transLengthDict,startCodonCoorDict,stopCodonCoorDict,options.left_length,options.right_length)
	plot_periodicity(start_density,stop_density,specific_counts,options.output_prefix,options.left_length,options.right_length)
	print("Finish the step of plot periodicity",file=sys.stderr)


if __name__=="__main__":
	parse_args()
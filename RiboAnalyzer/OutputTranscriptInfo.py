#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
@Author: Li Fajin
@since: 2019-08-12 09:57:47
@lastTime: 2019-08-12 15:49:18
@LastAuthor: Li Fajin
@Description:
	This script is used for extracting the longest and all transcript of each gene and output its
	annotation information.There are four inputs:
	1) coorFile: generated by RiboCode --> transcripts_cds.txt
	2) transcriptFile: generated by RiboCode --> transcripts_sequences.txt
	3) gtfFile: annotation file with gtf format for specific species
	4) TransInfoFile: the name of output file
	Usage: python OutputTranscriptInfo.py -c transcripts_cds.txt -g gtfFile.gtf -f transcripts_sequences.fa -o longest.transcripts.info.txt -O all.transcripts.info.txt
'''


from .FunctionDefinition import *
from itertools import groupby




def main():
    parser=create_parser_for_output_transInfo()
    (options,args)=parser.parse_args()
    if options.longestTransInfo and not options.allTransInfo:
        print("Starting outputing longest trans...")
        get_longest_transcripts_information(options.coorFile,options.transcriptFile,options.gtfFile,options.longestTransInfo)
        print("Finishing!")
    elif options.allTransInfo and not options.longestTransInfo:
        print("Starting outputing all trans...")
        get_all_transcripts_information(options.coorFile,options.transcriptFile,options.gtfFile,options.allTransInfo)
        print("Finishing!")
    elif options.longestTransInfo and options.allTransInfo:
        print("Starting outputing longest trans...")
        get_longest_transcripts_information(options.coorFile,options.transcriptFile,options.gtfFile,options.longestTransInfo)
        print("Finishing!")
        print("Starting outputing all trans...")
        get_all_transcripts_information(options.coorFile,options.transcriptFile,options.gtfFile,options.allTransInfo)
        print("Finishing!")
    else:
        raise IOError("Input error! check your -O or -o parameters.")


if __name__=="__main__":
	main()

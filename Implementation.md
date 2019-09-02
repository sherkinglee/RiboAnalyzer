# **RiboAnalyzer: a python package for downstream analysis of ribosome profiling data**
The **RiboAnalyzer** is a python package used for downstream analysis of ribosome profiling data. This package has four function parts:

+ **Quality Control (QC)**: Quality control for ribosome profiling data, containing periodicity checking, reads distribution among different reading frames,length distribution of ribosome footprints and DNA contaminations.
+ **Metagene Analysis (MA)**: Metagene analysis among different samples to find possible ribosome stalling events.
+ **Feature Analysis (FA)**: Feature analysis among different gene sets identified in MA step to explain the possible ribosome stalling.
+ **Enrichment Analysis (EA)**: Enrichment analysis to find possible co-translation events.

In this file, we will show you how to use our **RiboAnalyzer** based on some published datasets.

# **Data preparation**

## **Datasets download**

we will use **[GSE89704](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89704)** and **[GSE116570](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116570)** datasets for following analysis. Genome sequence file and annotation file were downloaded from **[ensemble](http://www.ensembl.org/Saccharomyces_cerevisiae/Info/Index)**.
```
## for GSE89704
for i in SRR50081{34..37};do fastq-dump $i;done

## for GSE116570
for i in SRR74712{49..56};do fastq-sump $i;done
```
## **Prepare annotation file**
+ **Prepare sequences and annotaiton files on transcriptome level.**
```
prepare_transcripts -g Saccharomyces.gtf -f Saccharomyces_genome.fa -o RiboCode_annot
```
Because there is almost no UTR regions for Saccharomyces cerevisiae, for better observation I extended 50 nt away from start codon and stop codon, which means coordinates in *gtf* file are increased about 50 nt for both start and end position.

+ **Prepare the longest transcript annotaion files.**
```
## command
OutputTranscriptInfo -c RiboCode_annot/transcripts_cds.txt -g Saccharomyces.gtf -f RiboCode_annot/transcripts_sequence.fa -o longest.transcripts.info.txt -O all.transcripts.info.txt
## output
Starting outputing longest trans...
6692  transcripts will be used in the follow analysis.

Finishing!
Starting outputing all trans...
```
+ **Prepare the sequence file for the longest transcripts**
```
GetProteinCodingSequence -i RiboCode_annot/transcripts_sequence.fa  -c longest.transcripts.info.txt -o <output_prefix> --mode whole --table 1
##
6692  transcripts will be used in the follow analysis.

Notes: There are 0 transcripts whose cds sequence cannot be divided by 3!
Finish the step of extracting sequences!
```
# **Preprocessing**
The **Preprocessing** contains some basic steps, such as quality control, adapters trimming, quality filtering, removing rRNA contamination, reads mapping, et al. Here we used part of samples from **[GSE89704](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89704)** as an example.
## **Quality control**
```
#!/bin/bash
#BSUB -J beforeQC.sh
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -n 1
#BSUB -q TEST-A
workdir=/Share2/home/lifj/Projects/03.RiboAnalyzer/2017_eIF5A/01.beforeQC
fastaFile=/Share2/home/lifj/Projects/03.RiboAnalyzer/2017_eIF5A/00.rawdata
for i in  SRR5008135  SRR5008136  SRR5008137  SRR5008134;
do
    fastqc $fastaFile/$i.fastq -o $workdir
done
```
## **Adapter trimming**
```
#!/bin/bash
#BSUB -J cutadaptRPF.sh
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -n 4
#BSUB -q TEST-A
workdir=/Share2/home/lifj/Projects/03.RiboAnalyzer/2017_eIF5A/02.cutadapt
fastaFile=/Share2/home/lifj/Projects/03.RiboAnalyzer/2017_eIF5A/00.rawdata
adapter=CTGTAGGCACCATCAAT
for i in  SRR5008135  SRR5008136  SRR5008137  SRR5008134;
do
        cutadapt -m 15 -M 35 --match-read-wildcards -a $adapter -o $workdir/$i.trimmed.fastq $fastaFile/$i.fastq > $workdir/$i_trimmed.log
done
```
## **Quality filtering**
```
#!/bin/bash
#BSUB -J filter.sh
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -n 1
#BSUB -q TEST-A
workdir=/Share2/home/lifj/Projects/03.RiboAnalyzer/2017_eIF5A/03.filter
fastaFile=/Share2/home/lifj/Projects/03.RiboAnalyzer/2017_eIF5A/02.cutadapt
for i in  SRR5008135  SRR5008136  SRR5008137  SRR5008134;
do
    fastq_quality_filter -Q33 -v -q 25 -p 75 -i $fastaFile/$i.trimmed.fastq -o $workdir/$i.trimmed.Qfilter.fastq > $workdir/$i.Qfilter.log
done
```
## **Remove rRNA contamination**
```
#!/bin/bash
bowtie_noncoding_index=/Share2/home/lifj/Reference/yeast/saccharomyces/Bowtie_noncoding_index/saccharomyces_noncoding
fastqFile=/Share2/home/lifj/Projects/03.RiboAnalyzer/2017_eIF5A/03.filter
workdir=/Share2/home/lifj/Projects/03.RiboAnalyzer/2017_eIF5A/05.contam
for i in  SRR5008135  SRR5008136  SRR5008137  SRR5008134;
do
    bsub -q TEST-A -o $i.out -e $i.err "bowtie -n 0 -y -a --norc --best --strata -S -p 4 -l 15 --un=$workdir/noncontam_$i.fastq $bowtie_noncoding_index -q $fastqFile/$i.trimmed.Qfilter.fastq $workdir/$i.alin"
done
```
## **Quality control**

This step is the same as the last one **Quality control** based on *fastq* files generated on the **Remove rRNA contamination** step.

## **Reads mapping**
```
## use one sample for example
#!/bin/bash
#BSUB -J STAR.sh
#BSUB -n 8
STAR_genome_index=/Share2/home/lifj/Reference/yeast/saccharomyces/STAR_genome_extended_index
workdir=/Share2/home/lifj/Projects/03.RiboAnalyzer/2017_eIF5A/07.STAR
fastqFile=/Share2/home/lifj/Projects/03.RiboAnalyzer/2017_eIF5A/05.contam
i=SRR5008134
mkdir -p $workdir/${i}_STAR

## mapping
STAR --runThreadN 8 --outFilterType Normal --outWigType wiggle --outWigStrand Stranded --outWigNorm RPM --alignEndsType EndToEnd --outFilterMismatchNmax 1 --outFilterMultimapNmax 1 --genomeDir $STAR_genome_index --readFilesIn $fastqFile/noncontam_$i.fastq --outFileNamePrefix  $workdir/${i}_STAR/$i. --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outSAMattributes All
## sort
samtools sort -T $workdir/${i}_STAR/$i.Aligned.toTranscriptome.out.sorted -o $workdir/${i}_STAR/$i.Aligned.toTranscriptome.out.sorted.bam $workdir/${i}_STAR/$i.Aligned.toTranscriptome.out.bam
## index
samtools index  $workdir/${i}_STAR/$i.Aligned.toTranscriptome.out.sorted.bam
```

# **RiboAnalyzer**
## **Quality Control (QC)**
+ **Checking 3-nt periodicity**
```
## use RiboCode:metaplots (recommeded)
#!/bin/bash
#BSUB -q TEST-A
workdir=/Share2/home/lifj/Projects/03.RiboAnalyzer/2017_eIF5A/08.periodicity
BamDir=/Share2/home/lifj/Projects/03.RiboAnalyzer/2017_eIF5A/07.STAR
RiboCode_annot=/Share2/home/lifj/Reference/yeast/saccharomyces/RiboCode/RiboCode_annotate_extend
for i in  SRR5008135  SRR5008136  SRR5008137  SRR5008134;
do
    metaplots -a $RiboCode_annot -r $BamDir/${i}_STAR/$i.Aligned.toTranscriptome.out.bam -o $workdir/$i
done

## use RiboAnalyzer
#!/bin/bash
#BSUB -J run_Periodicity.sh
#BSUB -q TEST-A
workdir=/Share2/home/lifj/Projects/03.RiboAnalyzer/2017_eIF5A/08.periodicity
bamFiles=/Share2/home/lifj/Projects/03.RiboAnalyzer/2017_eIF5A/07.STAR
Ref=/Share2/home/lifj/Reference/yeast/saccharomyces
for i in  SRR5008135  SRR5008136  SRR5008137  SRR5008134;do
    Periodicity -i $bamFiles/${i}_STAR/$i.Aligned.toTranscriptome.out.sorted.bam -o $workdir/$i -c $Ref/longest.transcripts.info.txt -L 25 -R 35
done
```
based on which construct **attributes.txt** file:
```
bamFiles    readLengths Offsets bamLegends
./SRR5008134.bam    27,28,29,30 11,12,13,14 si-Ctrl-1
./SRR5008135.bam    27,28,29,30 11,12,13,14 si-Ctrl-2
./SRR5008136.bam    27,28,29    11,12,13    si-eIF5A-1
./SRR5008137.bam    27,28,29    11,12,13    si-eIF5A-2
```
![RiboCode](F:\lifj\python_test\Ribo-seq\RiboAnalyzer\Implementation\periodicity\RiboCode.png)
![RiboAnalyzer](F:\lifj\python_test\Ribo-seq\RiboAnalyzer\Implementation\periodicity\RiboAnalyzer.png)
![RiboAnalyzer-total]()

+ **Reads distribution among different reading frames.**
```
RiboDensityOfDiffFrames -f attributes.txt -c longest.transcripts.info.txt -o <output_prefix>  --plot yes
```
![plot](#)


+ **Length distribution**
```
LengthDistribution -i SRR5008135.trimmed.Qfilter.fastq -o SRR5008135
```
![Length](#)

+ **Checking DNA contamination**.
```
StatisticReadsOnDNAsContam -i  SRR5008135.Aligned.sortedByCoord.out.bam  -g Saccharomyces.gtf  -o  <output_prefix>
```

![DNA](#)

## **Metagene Analysis**
+ **Metagene analysis along the whole transcript region.**
```
MetageneAnalysisForTheWholeRegions -f attributes.txt -c longest.transcripts.info.txt -o <output_prefix> -b 15,90,15 -l 100 -n 10 -m 1 -e 5
```
![plot](#)

```
PlotMetageneAnalysisForTheWholeRegions -i <output_prefix_scaled_density_dataframe.txt> -o <output_prefix> -g si-Ctrl,si-eIF5A -r si-Ctrl-1,si-Ctrl-2__si-eIF5A-1,si-eIF5A-2 -b 15,90,15 --mode all --xlabel-loc -0.4
```

+ **Metagene analysis on CDS regions.**
```
MetageneAnalysis -f attributes.txt -c longest.transcripts.info.txt -o <output_prefix> -U codon -M RPKM -u 0 -d 500 -l 100 -n 10 -m 1 -e 5 --norm yes -y 100 --CI 0.95 --type CDS
```

![plot](#)

```
PlotMetageneAnalysis -i <output_prefix_CDS_normed_dataframe.txt>  -o <output_prefix> -u 0 -d 500 -g si-Ctrl,si-eIF5A -r si-Ctrl-1,si-Ctrl-2__si-eIF5A-1,si-eIF5A-2 -U codon --CI 0.95
```

+ **Metagene analysis on UTR regions.**
```
MetageneAnalysis -f attributes.txt -c longest.transcripts.info.txt -o <output_prefix> -U nt -M RPKM -u 50 -d 50 -l 50 -n 10 -m 1 -e 5 --norm yes -y 100 --CI 0.95 --type UTR
```
![plot](#)
```
PlotMetageneAnalysis -i <output_prefix_UTR_dataframe.txt>  -o <output_prefix> -u 50 -d 50 -g si-Ctrl,si-eIF5A -r si-Ctrl-1,si-Ctrl-2__si-eIF5A-1,si-eIF5A-2 -U nt --CI 0.95
```

+ **Polarity calculation.**
```
PolarityCalculation -f attributes.txt -c longest.transcripts.info.txt -o <output_prefix> -n 64
```
![plot](#)

```
PlotPolarity -i <output_prefix_polarity_dataframe.txt> -o <output_prefix> -g si-Ctrl,si-eIF5A -r si-Ctrl-1,si-Ctrl-2__si-eIF5A-1,si-eIF5A-2  -y 5
```

## **Feature Analysis (FA)**
+ **Pick out transcripts enriched ribosomes on specific region.**
```
RiboDensityForSpecificRegion -f attributes.txt -c longest.transcripts.info.txt -o <output_prefix> -U codon -M RPKM -L 1 -R 100
```

+ **Ribosome density at each kind of AA or codon.**
    ```
    RiboDensityAtEachKindAAOrCodon -f <attributes.txt> -c <longest.transcripts.info.txt> -o <output_prefix> -M RPKM -S <select_trans.txt> -l 100 -n 10 --table 1 -F <longest_cds_sequence.fa>
    ```

    The output of this step is a file containing ribosome density at each kind of codon or amino acid like this:
    ```
    AA      codon   ctrl-1       ctrl-2       treat-1      treat-2
    K       AAA     0.026594288206047208    0.029364263341463245    0.015578462857160988    0.014649362490996175
    N       AAC     0.023361565164523757    0.025593616862523674    0.025286390202746506    0.025597799715113747
    ```
    We can also calculate the ribosome density of each kind of AA or codon at a specific region, for example codon 25-75 by add *-u* and *-d* parameters.
    ```
    RiboDensityAtEachKindAAOrCodon -f <attributes.txt> -c <longest.transcripts.info.txt> -o <output_prefix> -M RPKM -S <select_trans.txt> -l 100 -n 10 --table 1 -F <longest_cds_sequence.fa> -u 1 -d 100
    ```

+ **Ribosome density on amino acids with positive or negative charge**
    ```
    PlotRiboDensityAtEachKindAAOrCodon -i <output_prefix_all_codon_density.txt> -o <output_prefix> -g si-Ctrl,si-eIF5A -r si-Ctrl-1,si-Ctrl-2__si-eIF5A-1,si-eIF5A-2 --level AA
    ```
+ **Pausing score of each triplete amino acid**.
    ```
    PausingScore -f <attributes.txt> -c <longest.transcripts.info.txt> -o <output_prefix> -M RPKM -S <select_trans.txt>  -l 100 -n 10 --table 1 -F  <longest_cds_sequence.fa>
    ```

    ```
    ProcessPausingScore -i <group1-rep1_pausing_score.txt,group1-rep2_pausing_score.txt,group2-rep1_pausing_score.txt,group2-rep2_pausing_score.txt> -o <output_prefix> -g <group1,group2> -r <group1-rep1,group2-rep2__group2-rep1,group2-rep2> --mode [raw|ratio] --ratio_filter 2 --pausing_score_filter 10
    ```

    Where the *0,1,2* represent the E,P,A site of  a tri-AA motif. Using **[Seq2Logo](http://www.cbs.dtu.dk/biotools/Seq2Logo/)** for logo plot:

    ```
    ## this tool needs a python2.7 environment
    PATH=/tools/python2:$PATH
    Seq2Logo.py -f output_prefix_pwm.txt -u probability -I 5 -o <output_preifx> --format PDF
    ```
+ **Ribosome density around the triplete amino acid (tri-AA) motifs**.

    +   **As for a specific tri-AA motif, such as poly-proline (PPP)**
    ```
    RiboDensityAroundTripleteAAMotifs -f <attributes.txt> -c <longest.transcripts.info.txt> -o PPP -M RPKM -S <select_trans.txt> -l 100 -n 10 --table 1 -F <longest_cds_sequence.fa> --type2 PPP --type1 PP
    ```

    ```
    PlotRiboDensityAroundTriAAMotifs -i PPP_motifDensity_dataframe.txt -o PPP -g si-Ctrl,si-eIF5A -r si-Ctrl-1,si-Ctrl-2__si-eIF5A-1,si-eIF5A-2 --mode mean --ymax 0.2
    ```
    + **As for a tri-AA motifs list which contains part of information like this**:
    ```
    motifs
    PPP
    DDD
    SPP
    LPP
    LPP
    ```
    Using following command to do such job:

    ```
    RiboDensityAroundTripleteAAMotifs -f <attributes.txt> -c <longest.transcripts.info.txt> -o <output_prefix> -M RPKM -S <select_trans.txt> -l 100 -n 10 --table 1 -F <longest_cds_sequence.fa> --motifList1 tri_AA_motifs1.txt --motifList2 tri_AA_motifs2.txt
    ```

    ```
    PlotRiboDensityAroundTriAAMotifs -i tri_AA_motifDensity_dataframe.txt -o output_prefix -g si-Ctrl,si-eIF5A -r si-Ctrl-1,si-Ctrl-2__si-eIF5A-1,si-eIF5A-2 --mode mean --ymax 0.2
    ```

+ **RPFdist calculation**.
    ```
    RPFdist -f <attributes.txt> -c <longest.transcripts.info.txt> -o <output_prefix> -M RPKM -S <select_trans.txt> -l 100 -n 10 -m 1 -e 5
    ```

+ **Local tRNA adaptation index and tRNA adaptation index**
    ```
    tAI -i 2954_up_cds_sequences.fa,1598_unblocked_cds_sequences.fa,433_down_cds_sequences.fa -N tRNA_GCNs_Saccharomyces_cerevisiae.txt -o  test -u 0 -d 500 -t 2954_up,1598_unblocked,433_down
    ```

    ```
    tAIPlot -i <output_prefix_tAI_dataframe.txt> -o <output_prefix> -u 0 -d 500 --mode all --start 5 --window 7 --step 1
    ```

+ **Local codon adaptation index and codon adaptation index**
```
cAI -i $workdir/2954_up_cds_sequences.fa,$workdir/1598_unblocked_cds_sequences.fa,$workdir/433_down_cds_sequences.fa -o  $workdir/featureAnalysis/test.yeast -u 0 -d 500 -t 2954_up,1598_unblocked,433_down --reference reference.fa
```

```
tAIPlot -i <output_prefix_local_cAI_dataframe.txt> -o <output_prefix> -u 0 -d 500 --mode all --start 5 --window 7 --step 1
```
+ **Hydrophobicity calculation and Charge amino acids**
```
hydropathyCharge -i $workdir/2954_up_cds_sequences.fa,$workdir/1598_unblocked_cds_sequences.fa,$workdir/433_down_cds_sequences.fa -t 2954_up,1598_unblocked,433_down -o $workdir/featureAnalysis/test.yeast_hydropathy -u 0 -d 500 --index $workdir/featureAnalysis/hydropathy.txt

hydropathyCharge -i $workdir/2954_up_cds_sequences.fa,$workdir/1598_unblocked_cds_sequences.fa,$workdir/433_down_cds_sequences.fa -t 2954_up,1598_unblocked,433_down -o $workdir/featureAnalysis/test.yeast_charge -u 0 -d 500 --index $workdir/featureAnalysis/AA_charge.txt


```

    the both index file looks like:
    ```
    ## hydrophobicity index download from AAindex
    AA	amio_acids	hydropathy
    A	Ala 	1.8
    R	Arg 	-4.5
    N	Asn 	-3.5
    D	Asp 	-3.5
    C	Cys 	2.5
    ...
    ## charge index
    AA	amio_acids	hydropathy
    A	Ala 	0
    R	Arg 	1
    N	Asn 	0
    D	Asp 	-1
    C	Cys 	0
    Q	Gln 	0
    E	Glu 	-1
    G	Gly 	0
    H	His 	1
    ...
    ```
    among which the hydrophobicity index are downloaded from [AAindex](https://www.genome.jp/aaindex/). Later on, use *PlotHydropathyCharge.py* for plot:

```
PlotHydropathyCharge -i test.yeast_hydropathy_values_dataframe.txt -o zzzz_hydropathy -u 0 -d 500 --mode all --ylab "Average Hydrophobicity"

PlotHydropathyCharge -i test.yeast_charge_values_dataframe.txt -o zzzz_charge -u 0 -d 500 --mode all --ylab "Average Charge"
```

## **Enrichment Analysis (EA)**
The enrichment analysis is used for identifying potential co-translation events just like [Ayala Shiber, et al](https://www.nature.com/articles/s41586-018-0462-y) did. Therefore, there are two input files for this step. One is the total translatome file and the other is the IP translatome file. The enrichment analysis contains four steps:
+ **Step 1: Calculate ribosome density at each position for each transcript**.
    ```
    RiboDensityAtEachPosition -c <longest.transcripts.info.txt> -f <attributes.txt> -o <output_prefix>  -U codon
    ```
    where *attributes.txt* and *longest.transcripts.info.txt* are files we have talked above. This step would generate two files for each sample. One is the ribosome density file which would be used for following analysis. The Other is the coverage file, containing coverage of each transcript on CDS region.
+ **Step 2: Calculate mean ribosome density for different replicates**.
    ```
    enrichmentMeanDensity.py -i rep1-density.txt,rep2-density.txt,rep3-density.txt -o <output_prefix>
    ```
    *-i*  represents ribosome density files of different replicates generated by *RiboDensityAtEachPosition.py*, which should be separated by comma. This would generate a mean density file for each sample like *output_prefix_mean_density.txt*, used for enrichment analysis. If there is no replicates, just pass this step.

+ **Step 3: Enrichment analysis**.
    ```
    ## all transcripts
    EnrichmentAnalysis --ctrl <total-translatome.txt> --treat <IP-translatome.txt> -c <longest.transcripts.info.txt> -o <output_prefix> -U codon -M RPKM -l 150 -n 10 -m 1 -e 30 --CI 0.95 -u 0 -d 500

    ## specific transcripts
    EnrichmentAnalysis --ctrl <total-translatome.txt> --treat <IP-translatome.txt> -c <longest.transcripts.info.txt> -o <output_prefix> -U codon -M RPKM -l 150 -n 10 -m 1 -e 30 --CI 0.95 -u 0 -d 500 -S <select_trans.txt>
    ```

+ **Step 4: Plot the enrichment ratio**.
    ```
    PlotEnrichmentRatio -i <enrichment_dataframe.txt> -o <output_prefix> -u 0 -d 500 --unit codon --mode all {--slide-window y --axvline 1 --axhline 60 --label legend_label}
    ```

+ **Notes: if you want to see the enrichment ratio for a single transcript, the *EnrichmentAnalysisForSingleTrans.py* would be helpful**.
    ```
    EnrichmentAnalysisForSingleTrans -i <output_prefix_codon_ratio.txt> -s <transcript_name> -o <output_prefix> -c <longest.trans.info.txt>  --id-type transcript_id --slide-window y --axhline 1
    ```
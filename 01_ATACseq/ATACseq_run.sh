#!/usr/bin/bash


path_root=/path_to/ATACseq_data
data_dir=$1
SampleID=$2
refGenome_bw2=/path_to/Data_ref/Mus_musculus/iGenomes/mm10/Bowtie2Index/genome
refGenome=/path_to/Data_ref/Mus_musculus/iGenomes/mm10/WholeGenomeFasta/genome.fa

mkdir -p ${path_root}/02_ATACrun/{fastqc,logINFO,cutadapt,alignment,unAssign,CoverageBW,peakCalling}

## fastqc pre
fastqc ${path_root}/${data_dir}/${SampleID}/${SampleID}_1.fq.gz -o ${path_root}/02_ATACrun/fastqc
fastqc ${path_root}/${data_dir}/${SampleID}/${SampleID}_2.fq.gz -o ${path_root}/02_ATACrun/fastqc

## cut adapt
cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -a GGGGGGGGGGGGX -A GGGGGGGGGGGGX -m 18 -q 20,20 --max-n=0.05 -e 0.2 -n 2 -o ${path_root}/02_ATACrun/cutadapt/${SampleID}_cut_1.fq.gz -p ${path_root}/02_ATACrun/cutadapt/${SampleID}_cut_2.fq.gz ${path_root}/${data_dir}/${SampleID}/${SampleID}_1.fq.gz ${path_root}/${data_dir}/${SampleID}/${SampleID}_2.fq.gz > ${path_root}/02_ATACrun/logINFO/${SampleID}_cutINFO.log

## fastqc after
fastqc ${path_root}/02_ATACrun/cutadapt/${SampleID}_cut_1.fq.gz -o ${path_root}/02_ATACrun/fastqc
fastqc ${path_root}/02_ATACrun/cutadapt/${SampleID}_cut_2.fq.gz -o ${path_root}/02_ATACrun/fastqc

## alignment
bowtie2 -p 24 --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -x ${refGenome_bw2} -1 ${path_root}/02_ATACrun/cutadapt/${SampleID}_cut_1.fq.gz -2 ${path_root}/02_ATACrun/cutadapt/${SampleID}_cut_2.fq.gz -S ${path_root}/02_ATACrun/alignment/${SampleID}.sam --un-conc ${path_root}/02_ATACrun/unAssign/ 2> ${path_root}/02_ATACrun/logINFO/${SampleID}_bowtie2.log 

## add head INFO
picard SortSam I=${path_root}/02_ATACrun/alignment/${SampleID}.sam O=${path_root}/02_ATACrun/alignment/${SampleID}.sort.sam SORT_ORDER=coordinate
samtools addreplacerg -r "@RG\tID:RG1\tSM:$SampleID\tPL:Illumina\tLB:Library.fa" -o ${path_root}/02_ATACrun/alignment/${SampleID}.sort.ChangePG.sam ${path_root}/02_ATACrun/alignment/${SampleID}.sort.sam
# picard MarkDuplicates -I ${path_root}/02_ATACrun/alignment/${SampleID}.sort.ChangePG.sam -O ${path_root}/02_ATACrun/alignment/${SampleID}.markDup.bam -M ${path_root}/02_ATACrun/alignment/${SampleID}.markDup.txt
## remove duplicates
picard MarkDuplicates -I ${path_root}/02_ATACrun/alignment/${SampleID}.sort.ChangePG.sam -O ${path_root}/02_ATACrun/alignment/${SampleID}.rmDup.bam --REMOVE_DUPLICATES true -M ${path_root}/02_ATACrun/alignment/${SampleID}.rmDup.txt

## remove MT reads
samtools view -h ${path_root}/02_ATACrun/alignment/${SampleID}.rmDup.bam | grep -v 'chrM' | samtools view -bS -o ${path_root}/02_ATACrun/alignment/${SampleID}.rmMt.bam

## filter 
samtools view -@ 8 -h -F 1804 -q 30 -b ${path_root}/02_ATACrun/alignment/${SampleID}.rmMt.bam -o ${path_root}/02_ATACrun/alignment/${SampleID}.uniq.bam
## sort and index
samtools sort -@ 8 ${path_root}/02_ATACrun/alignment/${SampleID}.uniq.bam > ${path_root}/02_ATACrun/alignment/${SampleID}.sort.bam
samtools index -@ 8 ${path_root}/02_ATACrun/alignment/${SampleID}.sort.bam

## convert to bw
bamCoverage -b ${path_root}/02_ATACrun/alignment/${SampleID}.sort.bam -o ${path_root}/02_ATACrun/CoverageBW/${SampleID}.sort.bw --binSize 50 --normalizeUsing RPKM --extendReads

## Peak Calling
macs2 callpeak -t ${path_root}/02_ATACrun/alignment/${SampleID}.sort.bam -g mm -f BAMPE -n ${SampleID} -q 0.05 --keep-dup=all --outdir ${path_root}/02_ATACrun/peakCalling/ 

## FRiP cal
READS_IN_PEAKS=$(bedtools intersect -abam ${path_root}/02_ATACrun/alignment/${SampleID}.sort.bam -b ${path_root}/02_ATACrun/peakCalling/${SampleID}_peaks.narrowPeak -u | wc -l)
TOTAL_READS=$(samtools view -c ${path_root}/02_ATACrun/alignment/${SampleID}.sort.bam)
FRiP=$(awk "BEGIN {printf \"%.4f\", $READS_IN_PEAKS / $TOTAL_READS}")

echo $FRiP > ${path_root}/02_ATACrun/logINFO/${SampleID}_FRiP.log

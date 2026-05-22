# #!/usr/bin/bash

# ### pre-processing
# fastqc

# multiqc


cores=32
projPath="/path_to/CUTTag_data"
ref="/path_to/Data_ref/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome"
spikeInRef="/path_to/Data_ref/Ecoil_Bowtie2Index/genome"
chromSize="/path_to/Data_ref/Mus_musculus/UCSC/mm10/mm10.chrom.sizes"

picardCMD="picard"
minQualityScore=2
binLen=500

mkdir -p ${projPath}/03_CUTTAG/02_DataQC/trim_galoredata
mkdir -p ${projPath}/03_CUTTAG/03_alignment/sam/bowtie2_summary
mkdir -p ${projPath}/03_CUTTAG/03_alignment/bam
mkdir -p ${projPath}/03_CUTTAG/03_alignment/bed
mkdir -p ${projPath}/03_CUTTAG/03_alignment/bedgraph
mkdir -p ${projPath}/03_CUTTAG/03_alignment/removeDuplicate/picard_summary
mkdir -p ${projPath}/03_CUTTAG/03_alignment/sam/fragmentLen
mkdir -p ${projPath}/03_CUTTAG/04_peakCalling/SEACR
mkdir -p ${projPath}/03_CUTTAG/03_alignment/bedgraph

cat ${projPath}/03_CUTTAG/01_sampeINFO/Sample_list.tsv | while read id
do

    ### trim_galore
    trim_galore -q 25 --phred33 --stringency 3 --length 30  --paired ${projPath}/01.RawData/${id}/${id}_1.fq.gz ${projPath}/01.RawData/${id}/${id}_2.fq.gz --gzip -o ${projPath}/03_CUTTAG/02_DataQC/trim_galoredata

    ### Alignment
    ## bowtie2
    ## bowtie2-build path/to/hg38/fasta/hg38.fa /path/to/bowtie2Index/hg38
    bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${ref} -1 ${projPath}/03_CUTTAG/02_DataQC/trim_galoredata/${id}_1_val_1.fq.gz -2 ${projPath}/03_CUTTAG/02_DataQC/trim_galoredata/${id}_2_val_2.fq.gz -S ${projPath}/03_CUTTAG/03_alignment/sam/${id}_bowtie2.sam &> ${projPath}/03_CUTTAG/03_alignment/sam/bowtie2_summary/${id}_bowtie2.txt

    ### picard duplicate
    ## Sort by coordinate
    $picardCMD SortSam I=${projPath}/03_CUTTAG/03_alignment/sam/${id}_bowtie2.sam O=${projPath}/03_CUTTAG/03_alignment/sam/${id}_bowtie2.sorted.sam SORT_ORDER=coordinate

    ## mark duplicates
    samtools addreplacerg -r "@RG\tID:RG1\tSM:$id\tPL:Illumina\tLB:Library.fa" -o ${projPath}/03_CUTTAG/03_alignment/sam/${id}_bowtie2.sorted.ChangePG.sam ${projPath}/03_CUTTAG/03_alignment/sam/${id}_bowtie2.sorted.sam

    $picardCMD MarkDuplicates -I ${projPath}/03_CUTTAG/03_alignment/sam/${id}_bowtie2.sorted.ChangePG.sam -O ${projPath}/03_CUTTAG/03_alignment/removeDuplicate/${id}_bowtie2.sorted.dupMarked.sam -M ${projPath}/03_CUTTAG/03_alignment/removeDuplicate/picard_summary/${id}_picard.dupMark.txt

    ## remove duplicates
    $picardCMD MarkDuplicates -I ${projPath}/03_CUTTAG/03_alignment/sam/${id}_bowtie2.sorted.ChangePG.sam -O ${projPath}/03_CUTTAG/03_alignment/removeDuplicate/${id}_bowtie2.sorted.rmDup.sam --REMOVE_DUPLICATES true -M ${projPath}/03_CUTTAG/03_alignment/removeDuplicate/picard_summary/${id}_picard.rmDup.txt

    ## Extract the 9th column from the alignment sam file which is the fragment length
    samtools view -F 0x04 ${projPath}/03_CUTTAG/03_alignment/removeDuplicate/${id}_bowtie2.sorted.rmDup.sam  | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' >${projPath}/03_CUTTAG/03_alignment/sam/fragmentLen/${id}_fragmentLen.txt

    ### Alignment results filtering and file format conversion
    samtools view -q $minQualityScore ${projPath}/03_CUTTAG/03_alignment/removeDuplicate/${id}_bowtie2.sorted.rmDup.sam > ${projPath}/03_CUTTAG/03_alignment/sam/${id}_bowtie2.qualityScore${minQualityScore}.sam

    ## Filter and keep the mapped read pairs
    samtools view -bS -F 0x04 ${projPath}/03_CUTTAG/03_alignment/removeDuplicate/${id}_bowtie2.sorted.rmDup.sam > ${projPath}/03_CUTTAG/03_alignment/bam/${id}_bowtie2.mapped.bam

    ## Convert into bed file format
    bedtools bamtobed -i ${projPath}/03_CUTTAG/03_alignment/bam/${id}_bowtie2.mapped.bam -bedpe > ${projPath}/03_CUTTAG/03_alignment/bed/${id}_bowtie2.bed

    ## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
    awk '$1==$4 && $6-$2 < 1000 {print $0}' ${projPath}/03_CUTTAG/03_alignment/bed/${id}_bowtie2.bed >${projPath}/03_CUTTAG/03_alignment/bed/${id}_bowtie2.clean.bed

    ## Only extract the fragment related columns
    cut -f 1,2,6 ${projPath}/03_CUTTAG/03_alignment/bed/${id}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  >${projPath}/03_CUTTAG/03_alignment/bed/${id}_bowtie2.fragments.bed

    awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' ${projPath}/03_CUTTAG/03_alignment/bed/${id}_bowtie2.fragments.bed | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n  >${projPath}/03_CUTTAG/03_alignment/bed/${id}_bowtie2.fragmentsCount.bin${binLen}.bed


    bedtools genomecov -bg -i ${projPath}/03_CUTTAG/03_alignment/bed/${id}_bowtie2.fragments.bed -g $chromSize > ${projPath}/03_CUTTAG/03_alignment/bedgraph/${id}_bowtie2.fragments.bedgraph

    # ### Peak calling
    seacr="SEACR_1.3.sh"

    $seacr ${projPath}/03_CUTTAG/03_alignment/bedgraph/${id}_bowtie2.fragments.bedgraph 0.01 non stringent ${projPath}/03_CUTTAG/04_peakCalling/SEACR/${id}_seacr_top0.01.peaks

    $seacr ${projPath}/03_CUTTAG/03_alignment/bedgraph/${id}_bowtie2.fragments.bedgraph ${projPath}/03_CUTTAG/03_alignment/bedgraph/NC_bowtie2.fragments.bedgraph non stringent ${projPath}/03_CUTTAG/04_peakCalling/SEACR/${id}_seacr_control.peaks

done




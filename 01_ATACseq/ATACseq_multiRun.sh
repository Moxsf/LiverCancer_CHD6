#!/usr/bin/bash


tmp_fifofile="/tmp/$$.fifo"
mkfifo $tmp_fifofile
exec 6<>$tmp_fifofile

thread_num=3
for ((i=0;i<${thread_num};i++))
do
    echo >&6
done 


TSV_FILE=/path_to/ATACseq_data/SampleINFO.tsv

while IFS=$'\t' read -r Col1 Col2 Col3
do
    read -u6
    {
        DataDir=$Col1
        fileID=$Col2
        SampleID=$Col3
        bash /path_to/ATACseq_script/ATACseq_run.sh $DataDir $fileID
        echo >&6
    }&

done < $TSV_FILE

wait

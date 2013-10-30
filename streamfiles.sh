#!/bin/sh

FILELIST="$1"
PATTERNSFILE="$2"
INPUTMODE="$3"

while read line; 
do 
if [ "$line" = "" ] ; then
    echo "empty line.";
else
    FILENAME="/.mounts/labs/simpsonlab/data/rice_reads/fastq/long/seq_only/$line.fastq.dat"    
    qsub -cwd -b y -l h_vmem=4g -m beas -M mbarsky@oicr.on.ca ./count_one_file.sh $PATTERNSFILE "50" $FILENAME $INPUTMODE
fi
done < $FILELIST

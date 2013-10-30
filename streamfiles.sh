#!/bin/sh

FILELIST="$1"
PATTERNSFILE="$2"
INPUTMODE="$3"

while read line; 
do 
if [ "$line" = "" ] ; then
    echo "empty line.";
else
    FILENAME="/my_directory/$line.fasta"    
    qsub -cwd -b y -l h_vmem=4g ./count_one_file.sh $PATTERNSFILE "50" $FILENAME $INPUTMODE
fi
done < $FILELIST

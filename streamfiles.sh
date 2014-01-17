#!/bin/sh

FILELIST="$1"
PATTERNSFILE="$2"
INPUTMODE="$3"
K="$4"

while read line; 
do 
if [ "$line" = "" ] ; then
    echo "empty line.";
else
    FILENAME="/.mounts/labs/simpsonlab/users/mbarsky/hg/seq_only/$line"    
    qsub -cwd -b y -l h_vmem=8g ./count_one_file.sh $PATTERNSFILE $K $FILENAME $INPUTMODE
fi
done < $FILELIST

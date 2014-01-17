#!/bin/sh
PATTERNSFILE="$1"
K="$2"
INPUTMODE="$4"
FILETOPROCESS="$3"
./streamcountunzipped -p $PATTERNSFILE -r -k $K -m 4000 -i $INPUTMODE  -c $FILETOPROCESS

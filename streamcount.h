#ifndef STREAMCOUNT_H
#define STREAMCOUNT_H

#include "kmers_to_kwtree.h"
#include "count_kmers.h"

//this is the main of the program

//*****************
//Part 1. Preprocess set of k-mers into a keyword tree: 
//***************

int convertAllKmersIntoKWTreeReturnTree (FILE *kmersFP, int inputType, SC_INT k, int includeRC, SC_INT memoryMB, KWTreeBuildingManager *manager); //in-memory version

//**************************
//Part 2 - Counting k-mers in one input file
//***************************

int streamAndCountOneFile(KWTCounterManager *manager);

//************************
//Part 3. Adding up counts for k-mer and its reverse complement
//************************
int combineSubstringCountsIntoKmersCounts(SC_INT totalSubstringCounts,SC_INT *substringCounts,
    SC_INT totalKmers, KmerInfo *kmersInfo, 
    SC_INT *kmersCounts, int includeRC );

#endif

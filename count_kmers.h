#ifndef COUNT_KMERS_H
#define COUNT_KMERS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

#include <zlib.h>
#include "kseq.h"

#include "omp.h"

#include "common.h"
#include "keyword_tree.h"

//struct to hold buffered sequences and their length for multi-threading
typedef struct BufferCell{
	char *seq;
    int len;
}BufferCell;

//bookkeeping for performing count
typedef struct KWTCounterManager {
	SC_INT numberOfKWTreeNodes; //number of slots in the array of tree nodes
	SC_INT k; //k - length of each k-mer
	SC_INT numberOfKWTreeLeaves; //total unique patterns to search - we are going to store that many counters per file, counters start from 1 in this array
	KWTNode *KWTree; //keyword tree holding all patterns to be counted
	SC_INT *substringCounts; //holds resulting count of each pattern - meaningful value begins from 1 in this array
	FILE *inputFP; //the pointer to a current file where counting is performed
    int inputType; //type of an input - default 0 - FASTA, LINES (1), FILE (2)
    int numberOfThreads;
}KWTCounterManager;

#endif

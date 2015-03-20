#ifndef KMERS_TO_KWTREE_H
#define KMERS_TO_KWTREE_H

#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include "common.h"
#include "keyword_tree.h"

//produces an array of k-mers by reading provided file. This can be treated as each line as a separate input, or the entire file as one long string
int fillKmersArrayAndInfo(KWTreeBuildingManager* manager, FILE *inputFP, 
	char **kmers, KmerInfo *kmersInfo, SC_INT maxPossibleNumberOfKmers);

//collects stats about possible number of k-mers to count
int collectKmerInputStats(FILE *inputFP, SC_INT k, int inputType, SC_INT *estimatedNumberOfKmers);

//builds kw tree index from all k-mers
int convertAllKmersIntoKWTree (char *kmersFileName, int inputType, SC_INT k, int includeRC, SC_INT memoryMB);

//saves binary tree to file for future use
int  saveTreeAndMapping (char *kwtreeFileName, KWTNode *KWtree, SC_INT totalNodes, KmerInfo *kmersInfo, SC_INT totalKmers);
#endif

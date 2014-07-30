#ifndef KEYWORD_TREE_H
#define KEYWORD_TREE_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "common.h"

//bookkeeping for building KWtree for pattern set - including repeating k-mers and reverse complement of each k-mer, which is counted as the same k-mer
typedef struct KWTreeBuildingManager
{
	long maxSetSize;  //how many tree nodes can be held in a given memory
	SC_INT estimatedNumberOfLeaves; //estimated number of leaves - depends on include RC or not
	SC_INT estimatedNumberOfKWTreeNodes;	//how many nodes it needs at most - to allocate memory
	SC_INT actualNumberOfKWTreeNodes; //number of nodes in actual tree - known only after tree is built
	SC_INT actualNumberOfLeaves; //corresponds to number of unique substrings contained in the tree
	SC_INT maxNumberOfLeaves; //after reading actual patterns - know how many leaves can be in the tree
	SC_INT treeLeavesNum; //number of leaves in the tree eq. number of UNIQUE k-mers
	SC_INT k; //k - length of each k-mer
	SC_INT originalNumberOfKmers; //number of total k-mers, actual number of patterns is one less - since we start from 1
    SC_INT estimatedNumberOfKmers; //we need that in order to clean memory afterwards
	int inputType;			//0-lines, 1 -file
	int includeReverseComplement;	//whether to include reverse complement: default yes (1)
	KWTNode *KWtree; //keyword tree holding all patterns to be counted
	KmerInfo *kmersInfo; //string and mapping information about each k-mer
    char **kmers; //actual strings
}KWTreeBuildingManager;

//for breadth first traversal of the tree - to add suffix links
typedef struct Queue  
{
	SC_INT first;
	SC_INT freeSpot;
	SC_INT counter;
	SC_INT *nodePointers;
}Queue;

//construction: build a trie and add the suffix links
int buildKeywordTree (KWTreeBuildingManager *manager, char ** patterns, KmerInfo *patternsInfo);
int addSuffixLinks (KWTNode *tree, SC_INT totalNodes);

//streaming through the tree
int streamOneString(KWTNode* KWTree,char *input,int strlength,SC_INT *patternCounts);
int streamOneStringMT(KWTNode* KWTree,char *input,int strlength,SC_INT *patternCounts); //multi-threaded version
#endif

// Copyright 2014 OICR
// Written by Marina Barsky (mbarsky@oicr.on.ca)
// Released under the GPL

/**
This is an implementation of Aho-Corasick algorithm for searching for a set of patterns in the input text.
The algorithm is twisted for the problem of counting a set of k-mers in arbitrarily large input in one pass over the disk data.

The algorithm includes two main steps:
1. building a keyword tree (with failure links) from all input k-mers
2. using traversal of this tree according to input sequences to increment count for each k-mer
**/

#include "keyword_tree.h"

//Builds a keyword tree (digital trie) from all k-mers 
//also inserts into the tree a reverse complement for each k-mer
//this is performed in time linear in the total number of characters in all k-mers
int buildKeywordTree (KWTreeBuildingManager *manager, char ** patterns, KmerInfo *patternsInfo) {
	SC_INT i,c;
	SC_INT currentNodeID = 0;
	int currentCharAsINT;
	
	SC_INT *repetitionFirstOccurrence = (SC_INT *)calloc(manager->maxNumberOfLeaves+1,sizeof (SC_INT));
	char *rcPattern = (char *)calloc(manager->k+1,sizeof (char));
	SC_INT leavesCounter = 1;
	SC_INT treeSlotsNum = 1;

	for( i = 0; i < manager->originalNumberOfKmers; i++)	{
        currentNodeID = 0; //starting from the root
		c=0; //starting from the first char
		currentCharAsINT = getCharValue(patterns[i][c]);
		if(currentCharAsINT<0) {
			fprintf(stderr,"Why an invalid char in the input k-mers?\n");
			return EXIT_FAILURE;
		}

		//follow the path if exists
		while (manager->KWtree[currentNodeID].children[currentCharAsINT]>0 && c<manager->k)	{			
			currentNodeID = manager->KWtree[currentNodeID].children[currentCharAsINT];           
			c++; 
			if(c<manager->k)
				currentCharAsINT = getCharValue(patterns[i][c]);          
		}

		//there are 2 cases: either the entire pattern is already in the tree - we reached a leaf
		//in this case we mark it - duplicate: repeated = 1; 
		if(c == manager->k){
			patternsInfo[i].repeated=1;
			patternsInfo[i].counterID = - manager->KWtree[currentNodeID].children[0];

			//also need to mark as repeated the first occurrence
			patternsInfo[repetitionFirstOccurrence[- manager->KWtree[currentNodeID].children[0]]].repeated = 1;			
		}
		else {
			//add a new path for the remaining suffix
			for(;c < manager->k;c++)	{
				SC_INT nextSlotID = treeSlotsNum++;
                
				currentCharAsINT = getCharValue(patterns[i][c]);
				manager->KWtree[currentNodeID].children[currentCharAsINT] = nextSlotID;
				currentNodeID=nextSlotID;				
			}

			//avoid overflow
			if(leavesCounter >= manager->maxNumberOfLeaves)	{
				fprintf(stderr,"UNEXPECTED ERROR: number of leaves in the tree exceeded estimated number\n");
				return EXIT_FAILURE;
			}

            //mark a leaf node
			manager->KWtree[currentNodeID].children[0]=-leavesCounter; //pointer to a corresponding pattern
			patternsInfo[i].counterID = leavesCounter;

			//record first occurrence of this pattern, in case it repeats later
			repetitionFirstOccurrence[leavesCounter] = i ;			
			leavesCounter++;			
		}

		//now repeat the same insertion but with reverse complement - record that this is rc for the same pattern
		if(manager->includeReverseComplement) {
			if(produceReverseComplement(patterns[i],rcPattern)!=EXIT_SUCCESS)
				return EXIT_FAILURE;
			
            currentNodeID=0; //starting from the root
			c=0; //starting from the first char
			currentCharAsINT = getCharValue(rcPattern[c]);

			//follow the path if exists
			while (manager->KWtree[currentNodeID].children[currentCharAsINT]>0 && c<manager->k)	{			
				currentNodeID = manager->KWtree[currentNodeID].children[currentCharAsINT];                
				c++; 
				if(c < manager->k)
					currentCharAsINT = getCharValue(rcPattern[c]);
			}

			//there are 2 cases: either the entire pattern is already in the tree -
			//in this case its original pattern has been already marked as duplicate
			if(c == manager->k)	{
				patternsInfo[i].rcCounterID = - manager->KWtree[currentNodeID].children[0];
			}
			else {				
				//add a new path for the remaining suffix
				for(;c < manager->k;c++) {
					SC_INT nextSlotID = treeSlotsNum++;
                   
					currentCharAsINT = getCharValue(rcPattern[c]);
					manager->KWtree[currentNodeID].children[currentCharAsINT] = nextSlotID;
					currentNodeID=nextSlotID;				
				}

				if(leavesCounter >= manager->maxNumberOfLeaves)	{
					fprintf(stderr,"UNEXPECTED ERROR 2: number of leaves in the tree exceeded number of expected patterns with RC\n");
					return EXIT_FAILURE;
				}	

				//mark a leaf node
				manager->KWtree[currentNodeID].children[0]=-leavesCounter; //pointer to a corresponding pattern
				patternsInfo[i].rcCounterID = leavesCounter;			
			
				leavesCounter++;							
			}
		}        
	}	
	
	//set totalUniquePatterns - leavesCounter
	manager->actualNumberOfLeaves = leavesCounter;
	manager->actualNumberOfKWTreeNodes = treeSlotsNum;

	free(repetitionFirstOccurrence);
	free(rcPattern);
	
    //add suffix links for fast search
	if(addSuffixLinks (&(manager->KWtree[0]), manager->actualNumberOfKWTreeNodes)!=EXIT_SUCCESS)
		return EXIT_FAILURE;

	if(PRINT_KWTREE)	printKeywordTree(manager->KWtree, 0,0);
	
	return EXIT_SUCCESS;
}


//(BFT) - breadth first traversal:
//we take first element from the queue - it becomes parent
//we find all its children and add suffix links for them - to the parent's child
//we push children into priority queue to be processed later
int addSuffixLinks (KWTNode *tree, int totalNodes) { 
	int i;
	Queue queue;
    
    //allocate queue for all possible nodes
	if(!(queue.nodePointers  = (SC_INT *)calloc(totalNodes,sizeof (SC_INT))))	{
		fprintf(stderr,"Unable to allocate memory for queue of nodes to add suffix links to the keyword tree.\n");
		return EXIT_FAILURE;
	}
	queue.first = 0;
	queue.freeSpot =1;
	queue.counter = 1;

	tree[0].suffixLinkID = -1;  //there is no suffix link from the root node
	queue.nodePointers[0]  = 0; //push the root node into a queue - first points to the root node
	
	while(queue.counter > 0)	{
		KWTNode* popParent=&tree[queue.nodePointers[queue.first]];
		queue.first++;
		if(queue.first == totalNodes)		{
			fprintf(stderr,"Unexpected error in priority queue 1 - number of elements exceeded allocated space %d: queue.first=%d \n",totalNodes,totalNodes);
			free(queue.nodePointers);
			return EXIT_FAILURE;
		}
		queue.counter--;

		for( i = 0; i < SIGMA; i++)	{
			if(popParent->children[i] > 0) {
				SC_INT suffixLinkID = popParent->suffixLinkID;
				int sl_found=0;
				while(suffixLinkID !=-1 && !sl_found) {
					KWTNode *currentNode = &tree [suffixLinkID];
					if(currentNode->children[i] > 0) {
						sl_found=1;
						tree[popParent->children[i]].suffixLinkID = currentNode->children[i];
					}
					else
						suffixLinkID = currentNode->suffixLinkID;
				}

				if(!sl_found) { //no proper suffix of the current keyword is in the tree
					tree[popParent->children[i]].suffixLinkID = 0; //set suffix link to the root
				}

				//push the child into the queue - if it is not a leaf node
				if( !(tree[popParent->children[i]].children[0] < 0) ) { //not a leaf node				
					queue.nodePointers[queue.freeSpot] = popParent->children[i];

					//advance freespot
					queue.freeSpot++;

					if(queue.freeSpot >= totalNodes){
						fprintf(stderr,"Unexpected error in priority queue 2- number of elements exceeded allocated space %d\n",totalNodes);
						free(queue.nodePointers);
						return EXIT_FAILURE;
					}
					queue.counter++;
				}
			}
		}
	}
	
	free(queue.nodePointers);
	return EXIT_SUCCESS;
}

//this increments a count of the corresponding leaf (substring found in the kw tree)
//but here it checks if the value is not beeing updated by another thread
//ensures atomicity of write
void incrementKmerCount (SC_INT *patternCounts, SC_INT patternID)  {//multi-threaded version
    SC_INT old = patternCounts[patternID];
    SC_INT new = old+1;
    if(new >= MAX_COUNT)
        return;
    
    while(!__sync_bool_compare_and_swap(&patternCounts[patternID], old, new))    {  
        old = patternCounts[patternID];
        new = old+1;      
    }
}

//search for all substrings of a current input line in the tree by one pass through the input line
//and simultaneous traversal of the kw tree
//suffix links ensure that there are no more than 2 node traversal operations per each character of the input

//this is a multi-threaded version, which uses the above function to update counts
int streamOneStringMT(KWTNode* KWTree,char *input,int strlength,SC_INT *patternCounts) {
	SC_INT currentNodeID = 0; //start from the root
	SC_INT currentPositionInInput = 0; //start from the first character
	int found=0;
	SC_INT suffixLinkID;
	int currentChar;
	int invalidChar = 0;

	while (currentPositionInInput < strlength && !invalidChar)
	{
		currentChar = getCharValue(input[currentPositionInInput]);
		if(currentChar<0)
		{
            currentNodeID = 0; //re-start from the root
            currentPositionInInput++;
        }
		else
		{
			//case 1: there is a child currentChar out of a current node - we follow the path down the tree
			if(KWTree[currentNodeID].children[currentChar]>0)
			{
				currentNodeID=KWTree[currentNodeID].children[currentChar];
				if(KWTree[currentNodeID].children[0]<0) //leaf node - stores a negated pattern ID -  update counter
				{
					incrementKmerCount (patternCounts,-KWTree[currentNodeID].children[0]);                    
				}
				currentPositionInInput++;
			}
			else //case 2. no child currentChar out of a current node
			{
				suffixLinkID =KWTree[currentNodeID].suffixLinkID;
				found =0;
				while(suffixLinkID!=-1 && !found) //follows suffix links until finds outgoing edge for currentChar or reached the root - and there is no outgoing edge for currentChar
				{
					currentNodeID = suffixLinkID;
					if(KWTree[currentNodeID].children[currentChar]>0)	{
						currentNodeID=KWTree[currentNodeID].children[currentChar];
						if(KWTree[currentNodeID].children[0]<0)  {//leaf node - stores a negated pattern ID -  update counter
							incrementKmerCount (patternCounts,-KWTree[currentNodeID].children[0]);  
						}
						currentPositionInInput++;
						found=1;
					}
					else {
						suffixLinkID= KWTree[currentNodeID].suffixLinkID;  //follow up
					}
				}
			
				if(suffixLinkID == -1) { //reached the root and there is no appropriate child from the root - that means we need to start from the root and start from the next character
					currentNodeID=0;
					currentPositionInInput++;
				}			
			}
		}		
	}

	return EXIT_SUCCESS;
}

//this is the same function for a single-threaded execution
int streamOneString(KWTNode* KWTree,char *input,int strlength,SC_INT *patternCounts) {
	SC_INT currentNodeID = 0; //start from the root
	SC_INT currentPositionInInput = 0; //start from the first character
	int found=0;
	SC_INT suffixLinkID;
	int currentChar;
	int invalidChar = 0;

	while (currentPositionInInput < strlength && !invalidChar) {
		currentChar = getCharValue(input[currentPositionInInput]);
		if(currentChar<0) {
            currentNodeID = 0; //re-start from the root
            currentPositionInInput++;
        }
		else {
			//case 1: there is a child currentChar out of a current node - we follow the path down the tree
			if(KWTree[currentNodeID].children[currentChar]>0) {
				currentNodeID=KWTree[currentNodeID].children[currentChar];
				if(KWTree[currentNodeID].children[0]<0) {//leaf node - stores a negated pattern ID -  update counter				
					if(patternCounts[-KWTree[currentNodeID].children[0]]<MAX_COUNT)
						patternCounts[-KWTree[currentNodeID].children[0]]++;                    
				}
				currentPositionInInput++;
			}
			else { //case 2. no child currentChar out of a current node			
				suffixLinkID =KWTree[currentNodeID].suffixLinkID;
				found =0;
				while(suffixLinkID!=-1 && !found) {//follows suffix links until finds outgoing edge for currentChar or reached the root - and there is no outgoing edge for currentChar
					currentNodeID = suffixLinkID;
					if(KWTree[currentNodeID].children[currentChar]>0) {
						currentNodeID=KWTree[currentNodeID].children[currentChar];
						if(KWTree[currentNodeID].children[0]<0) {//leaf node - stores a negated pattern ID -  update counter						
							if(patternCounts[-KWTree[currentNodeID].children[0]]<MAX_COUNT)
								patternCounts[-KWTree[currentNodeID].children[0]]++;
						}
						currentPositionInInput++;
						found=1;
					}
					else {
						suffixLinkID= KWTree[currentNodeID].suffixLinkID;  //follow up
					}
				}
			
				if(suffixLinkID == -1) { //reached the root and there is no appropriate child from the root - that means we need to start from the root and start from the next character
					currentNodeID=0;
					currentPositionInInput++;
				}			
			}
		}		
	}

	return EXIT_SUCCESS;
}



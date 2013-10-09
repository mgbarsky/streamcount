#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "KeywordTree.h"

int buildKeywordTree (KWTreeBuildingManager *manager, char ** patterns, int64_t *totalPatterns, int64_t *totalUniquePatterns)
{
	int i,c;
	int currentNodeID = 0;
	
	INT currentCharAsINT;
	INT newPatternCounter;
	int numPatterns = *totalPatterns;
	char *currentPattern = (char *)calloc(manager->k+1,sizeof (char *));

	for(i=0;i<*totalPatterns;i++)
	{
		currentNodeID=0; //starting from the root
		c=0; //starting from the first char
		currentCharAsINT = getCharValue(patterns[i][c]);
		//follow the path if exists
		while (manager->KWtree[currentNodeID].children[currentCharAsINT]>0 && c<manager->k)
		{			
			currentNodeID = manager->KWtree[currentNodeID].children[currentCharAsINT];
			c++; 
			if(c<manager->k)
				currentCharAsINT = getCharValue(patterns[i][c]);
		}
		//there are 2 cases: either the entire pattern is already in the tree -
		//in this case we remove it - duplicate - by marking first character of a pattern -1
		if(c==manager->k)
		{
			patterns[i][0]=-1;
			numPatterns--;  //that is to test after removing duplicates
		}
		else
		{
			//add a new path for the remaining suffix
			for(;c < manager->k;c++)
			{
				INT nextSlotID = manager->treeSlotsNum++;
				currentCharAsINT = getCharValue(patterns[i][c]);
				manager->KWtree[currentNodeID].children[currentCharAsINT] = nextSlotID;
				currentNodeID=nextSlotID;
				
			}
			//mark a leaf node
			manager->KWtree[currentNodeID].children[0]=-(i+1); //pointer to a corresponding pattern
		}
	}
	
	newPatternCounter=1; //we will store patterns starting from 1 - to make use of negative numbers
	//traverse the tree and add a unique set of patterns - we will store them back into patterns array	
	if(traverseAndRecordPatterns(manager, currentPattern,0, patterns, &newPatternCounter,0)!=0)
		return 1;
	
	//sanity check
	if(newPatternCounter!=numPatterns+1)
	{
		printf("Unexpected error: number %d of unique patterns is not equal to expected number %d\n",newPatternCounter,numPatterns);
		return 1;
	}

	//set totalUniquePatterns - newPatternCounter
	*totalUniquePatterns=newPatternCounter;
	if(DEBUG_KWTREE)
	{
		printPatterns(patterns, newPatternCounter);
	}
	free(currentPattern);

	 
	if(addSuffixLinks (&(manager->KWtree[0]), manager->treeSlotsNum)!=0)
		return 1;

	if(DEBUG_KWTREE)
	{
		printKeywordTree(manager->KWtree, 0,0);
	}
	return 0;
}

//this takes one string of chars converted into INTs of alphabet, and streams this string through keyword tree
//when leaf node is reached - the count is incremented
//The streaming takes at most 2N operations, where N is the length of the input string
int streamOneString(KWTCounterManager *manager, INT *input, int length)
{
	int currentNodeID = 0; //start from the root
	int currentPositionInInput = 0; //start from the first character
	int found=0;
	int suffixLinkID;
	INT currentChar;

	while (currentPositionInInput < length)
	{
		currentChar = input[currentPositionInInput];
		//case 1: there is a child currentChar out of a current node - we follow the path down the tree
		if(manager->KWTree[currentNodeID].children[currentChar]>0)
		{
			currentNodeID=manager->KWTree[currentNodeID].children[currentChar];
			if(manager->KWTree[currentNodeID].children[0]<0) //leaf node - stores a negated pattern ID -  update counter
			{
				if(manager->patternCounts[-manager->KWTree[currentNodeID].children[0]]<MAX_COUNT)
					manager->patternCounts[-manager->KWTree[currentNodeID].children[0]]++;
			}
			currentPositionInInput++;
		}
		else //case 2. no child currentChar out of a current node
		{
			suffixLinkID =manager->KWTree[currentNodeID].suffixLinkID;
			found =0;
			while(suffixLinkID!=-1 && !found) //follows suffix link until finds outgoing edge for currentChar or reached the root - and there is no outgoing edge for currentChar
			{
				currentNodeID = suffixLinkID;
				if(manager->KWTree[currentNodeID].children[currentChar]>0)
				{
					currentNodeID=manager->KWTree[currentNodeID].children[currentChar];
					if(manager->KWTree[currentNodeID].children[0]<0) //leaf node - stores a negated pattern ID -  update counter
					{
						if(manager->patternCounts[-manager->KWTree[currentNodeID].children[0]]<MAX_COUNT)
							manager->patternCounts[-manager->KWTree[currentNodeID].children[0]]++;
					}
					currentPositionInInput++;
					found=1;
				}
				else
				{
					suffixLinkID= manager->KWTree[currentNodeID].suffixLinkID;  //follow up
				}
			}
			
			if(suffixLinkID == -1)  //reached the root and there is no appropriate child from the root - that means we need to start from the root and start from the next character
			{
				currentNodeID=0;
				currentPositionInInput++;
			}			
		}		
	}
	return 0;
}

//this recursively traverses (DFT) the keyword tree to collect all unique patterns, which are stored back into patterns array starting from position 1
//each leaf of the tree now stores a negated position of a pattern in this array at child[0]
int traverseAndRecordPatterns(KWTreeBuildingManager *manager, char *currentPattern, int posInPattern,
	char **patterns, INT *newPatternCounter, int parentNodeID)
{
	INT i;
	
	//we reached the leaf
	if(manager->KWtree[parentNodeID].children[0]<0)
	{
		//copy pattern to patterns array
		memcpy(patterns[*newPatternCounter],currentPattern,manager->k);
		patterns[*newPatternCounter][manager->k]='\0';

		//update leaf node id
		manager->KWtree[parentNodeID].children[0]=-(*newPatternCounter);
		(*newPatternCounter)++;
		
	}

	else
	{
		//loop through existing children
		for(i=0;i<SIGMA;i++)
		{
			if(manager->KWtree[parentNodeID].children[i] >0)
			{
				char currentChar = getCharFromINT(i);
				currentPattern[posInPattern] = currentChar;

				//DFS - to reach the leaf
				traverseAndRecordPatterns(manager, currentPattern, posInPattern+1,
					patterns, newPatternCounter, manager->KWtree[parentNodeID].children[i] );
			}
		}
	}
	return 0;
}


//(BFT) - breadth first traversal - 
//we take first element from the queue - it becomes parent
//we find all its children and add suffix links for them - to the parent's child
//we push children into priority queue to be processed later
int addSuffixLinks (KWTNode *tree, int totalNodes)
{ 
	int i;
	Queue queue;
	if(!(queue.nodePointers  =(INT *)calloc(totalNodes,sizeof (INT))))
	{
		printf("Unable to allocate memory for queue of nodes to add suffix links to the keyword tree.\n");
		return 1;
	}
	queue.first = 0;
	queue.freeSpot =1;
	queue.counter = 1;

	tree[0].suffixLinkID = -1;  //there is no suffix link from the root node
	queue.nodePointers[0]  = 0; //push the root node into a queue - first points to the root node
	
	while(queue.counter>0)
	{
		KWTNode* popParent=&tree[queue.nodePointers[queue.first]];
		queue.first++;
		if(queue.first==totalNodes)
		{
			printf("Unexpected error in priority queue - number of elements exceeded allocated space\n");
			free(queue.nodePointers);
			return 1;
		}
		queue.counter--;

		for(i=0;i<SIGMA;i++)
		{
			if(popParent->children[i]>0)
			{
				INT suffixLinkID = popParent->suffixLinkID;
				int sl_found=0;
				while(suffixLinkID !=-1 && !sl_found)
				{
					KWTNode *currentNode = &tree [suffixLinkID];
					if(currentNode->children[i]>0)
					{
						sl_found=1;
						tree[popParent->children[i]].suffixLinkID = currentNode->children[i];
					}
					else
						suffixLinkID = currentNode->suffixLinkID;
				}

				if(!sl_found) //no proper suffix of the current keyword is in the tree
				{
					tree[popParent->children[i]].suffixLinkID=0; //set suffix link to the root
				}

				//push the child into the queue - if it is not a leaf node
				if(!(tree[popParent->children[i]].children[0]<0)) //not a leaf node
				{
					queue.nodePointers[queue.freeSpot] = popParent->children[i];

					//advance freespot
					queue.freeSpot++;

					if(queue.freeSpot==totalNodes)
					{
						printf("Unexpected error in priority queue - number of elements exceeded allocated space\n");
						free(queue.nodePointers);
						return 1;
					}
					queue.counter++;
				}
			}
		}
	}
	
	free(queue.nodePointers);
	return 0;
}


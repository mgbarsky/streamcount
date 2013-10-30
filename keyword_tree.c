#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "KeywordTree.h"

int buildKeywordTree (KWTreeBuildingManager *manager, char ** patterns, KmerInfo *patternsInfo, int64_t *totalPatterns, int64_t *totalUniquePatterns)
{
	int i,c;
	int currentNodeID = 0;
	int maxNumberOfleaves = (manager->includeReverseComplement)? (*totalPatterns)*2:(*totalPatterns);
	maxNumberOfleaves+=2;  //that's because we start counting from 1

	INT currentCharAsINT;
	
	int numUniquePatterns = maxNumberOfleaves;
	int *repetitionFirstOccurrence = (int *)calloc(maxNumberOfleaves,sizeof (int));
	char *rcPattern = (char *)calloc(manager->k+1,sizeof (char));
	printf("Max number of leaves = %d\n",maxNumberOfleaves);

	for(i=0;i<*totalPatterns;i++)
	{
		currentNodeID=0; //starting from the root
		c=0; //starting from the first char
		currentCharAsINT = getCharValue(patterns[i][c]);
		if(currentCharAsINT<0)
		{
			printf("Why invalid char in the pattern?\n");
			return 1;
		}
		//follow the path if exists
		while (manager->KWtree[currentNodeID].children[currentCharAsINT]>0 && c<manager->k)
		{			
			currentNodeID = manager->KWtree[currentNodeID].children[currentCharAsINT];
			c++; 
			if(c<manager->k)
				currentCharAsINT = getCharValue(patterns[i][c]);
		}
		//there are 2 cases: either the entire pattern is already in the tree -
		//in this case we mark it - duplicate: repeated = 1; 
		if(c==manager->k)
		{
			patternsInfo[i].repeated=1;
			patternsInfo[i].counterID = - manager->KWtree[currentNodeID].children[0];

			//also need to mark as repeated the first occurrence
			patternsInfo[repetitionFirstOccurrence[- manager->KWtree[currentNodeID].children[0]]].repeated = 1;
			numUniquePatterns--;  //removing duplicates
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
			if(manager->treeLeavesNum >= maxNumberOfleaves)
			{
				printf("UNEXPECTED ERROR: number of leaves in the tree exceeded number of inserted patterns\n");
				return 1;
			}
			manager->KWtree[currentNodeID].children[0]=-manager->treeLeavesNum; //pointer to a corresponding pattern
			patternsInfo[i].counterID = manager->treeLeavesNum;
			//record first occurrence of this pattern, in case it repeats later
			repetitionFirstOccurrence[manager->treeLeavesNum] = i ;
			
			manager->treeLeavesNum++;
			
		}

		//now repeat the same insertion but with reverse complement
		if(manager->includeReverseComplement)
		{
			if(produceReverseComplement(patterns[i],rcPattern)!=0)
				return 1;
			currentNodeID=0; //starting from the root
			c=0; //starting from the first char
			currentCharAsINT = getCharValue(rcPattern[c]);
			//follow the path if exists
			while (manager->KWtree[currentNodeID].children[currentCharAsINT]>0 && c<manager->k)
			{			
				currentNodeID = manager->KWtree[currentNodeID].children[currentCharAsINT];
				c++; 
				if(c<manager->k)
					currentCharAsINT = getCharValue(rcPattern[c]);
			}
			//there are 2 cases: either the entire pattern is already in the tree -
			//in this case its original pattern has been already marked as duplicate
			if(c == manager->k)
			{
				patternsInfo[i].rcCounterID = - manager->KWtree[currentNodeID].children[0];
			}
			else
			{				
				//add a new path for the remaining suffix
				for(;c < manager->k;c++)
				{
					INT nextSlotID = manager->treeSlotsNum++;
					currentCharAsINT = getCharValue(rcPattern[c]);
					manager->KWtree[currentNodeID].children[currentCharAsINT] = nextSlotID;
					currentNodeID=nextSlotID;				
				}

				if(manager->treeLeavesNum >= maxNumberOfleaves)
				{
					printf("UNEXPECTED ERROR 2: number of leaves in the tree exceeded number of inserted patterns\n");
					return 1;
				}	

				//mark a leaf node
				manager->KWtree[currentNodeID].children[0]=-manager->treeLeavesNum; //pointer to a corresponding pattern
				patternsInfo[i].rcCounterID = manager->treeLeavesNum;			
			
				manager->treeLeavesNum++;
							
			}
		}
	}
	
	
	//set totalUniquePatterns - numUniquePatterns
	*totalUniquePatterns=numUniquePatterns;
	
	free(repetitionFirstOccurrence);
	free(rcPattern);
	 
	if(addSuffixLinks (&(manager->KWtree[0]), manager->treeSlotsNum)!=0)
		return 1;

	if(DEBUG_KWTREE)
	{
		printKeywordTree(manager->KWtree, 0,0);
	}
	return 0;
}

//this takes one string of chars, and streams this string through keyword tree
//when leaf node is reached - the count is incremented
//The streaming takes at most 2N operations, where N is the length of the input string
int streamOneStringUnchanged(KWTCounterManager *manager, char *input, int strlength)
{
	int currentNodeID = 0; //start from the root
	int currentPositionInInput = 0; //start from the first character
	int found=0;
	int suffixLinkID;
	int currentChar;
	int invalidChar = 0;
	while (currentPositionInInput < strlength && !invalidChar)
	{
		currentChar = getCharValue(input[currentPositionInInput]);
		if(currentChar<0)
			invalidChar=1;
		else
		{
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


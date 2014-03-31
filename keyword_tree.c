#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "General.h"
#include "StreamCount.h"

int buildKeywordTree (KWTreeBuildingManager *manager, char ** patterns, KmerInfo *patternsInfo)
{
	INT i,c;
	INT currentNodeID = 0;
	int currentCharAsINT;
	
	INT *repetitionFirstOccurrence = (INT *)calloc(manager->maxNumberOfLeaves+1,sizeof (INT));
	char *rcPattern = (char *)calloc(manager->k+1,sizeof (char));
	INT leavesCounter = 1;
	INT treeSlotsNum = 1;

    //fprintf(stderr,"Max number of leaves = %ld\n",manager->maxNumberOfLeaves);

	for(i=0;i<manager->originalNumberOfKmers;i++)
	{
        currentNodeID=0; //starting from the root
		c=0; //starting from the first char
		currentCharAsINT = getCharValue(patterns[i][c]);
		if(currentCharAsINT<0)
		{
			fprintf(stderr,"Why invalid char in the input k-mers?\n");
			return EXIT_FAILURE;
		}

        
		//follow the path if exists
		while (manager->KWtree[currentNodeID].children[currentCharAsINT]>0 && c<manager->k)
		{			
			currentNodeID = manager->KWtree[currentNodeID].children[currentCharAsINT];
           
			c++; 
			if(c<manager->k)
				currentCharAsINT = getCharValue(patterns[i][c]);
          
		}

		//there are 2 cases: either the entire pattern is already in the tree - we reached a leaf
		//in this case we mark it - duplicate: repeated = 1; 
		if(c==manager->k)
		{
			patternsInfo[i].repeated=1;
			patternsInfo[i].counterID = - manager->KWtree[currentNodeID].children[0];

			//also need to mark as repeated the first occurrence
			patternsInfo[repetitionFirstOccurrence[- manager->KWtree[currentNodeID].children[0]]].repeated = 1;			
		}
		else
		{
			//add a new path for the remaining suffix
			for(;c < manager->k;c++)
			{
				INT nextSlotID = treeSlotsNum++;
                
				currentCharAsINT = getCharValue(patterns[i][c]);
				manager->KWtree[currentNodeID].children[currentCharAsINT] = nextSlotID;
				currentNodeID=nextSlotID;				
			}

			//mark a leaf node
			if(leavesCounter >= manager->maxNumberOfLeaves)
			{
				fprintf(stderr,"UNEXPECTED ERROR: number of leaves in the tree exceeded estimated number\n");
				return EXIT_FAILURE;
			}

			manager->KWtree[currentNodeID].children[0]=-leavesCounter; //pointer to a corresponding pattern
			patternsInfo[i].counterID = leavesCounter;
			//record first occurrence of this pattern, in case it repeats later
			repetitionFirstOccurrence[leavesCounter] = i ;			
			leavesCounter++;			
		}

		//now repeat the same insertion but with reverse complement
		if(manager->includeReverseComplement)
		{
			if(produceReverseComplement(patterns[i],rcPattern)!=EXIT_SUCCESS)
				return EXIT_FAILURE;
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
					INT nextSlotID = treeSlotsNum++;
                   
					currentCharAsINT = getCharValue(rcPattern[c]);
					manager->KWtree[currentNodeID].children[currentCharAsINT] = nextSlotID;
					currentNodeID=nextSlotID;				
				}

				if(leavesCounter >= manager->maxNumberOfLeaves)
				{
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
	 
	if(addSuffixLinks (&(manager->KWtree[0]), manager->actualNumberOfKWTreeNodes)!=EXIT_SUCCESS)
		return EXIT_FAILURE;

	if(PRINT_KWTREE)
	{
		printKeywordTree(manager->KWtree, 0,0);
	}
	return EXIT_SUCCESS;
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
		fprintf(stderr,"Unable to allocate memory for queue of nodes to add suffix links to the keyword tree.\n");
		return EXIT_FAILURE;
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
			fprintf(stderr,"Unexpected error in priority queue 1 - number of elements exceeded allocated space %d: queue.first=%d \n",totalNodes,totalNodes);
			free(queue.nodePointers);
			return EXIT_FAILURE;
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

					if(queue.freeSpot>=totalNodes)
					{
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

int streamOneString(KWTNode* KWTree,char *input,int strlength,INT *patternCounts)
{
	INT currentNodeID = 0; //start from the root
	INT currentPositionInInput = 0; //start from the first character
	int found=0;
	INT suffixLinkID;
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
					if(patternCounts[-KWTree[currentNodeID].children[0]]<MAX_COUNT)
                    {
						patternCounts[-KWTree[currentNodeID].children[0]]++;                        
                    }
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
					if(KWTree[currentNodeID].children[currentChar]>0)
					{
						currentNodeID=KWTree[currentNodeID].children[currentChar];
						if(KWTree[currentNodeID].children[0]<0) //leaf node - stores a negated pattern ID -  update counter
						{
							if(patternCounts[-KWTree[currentNodeID].children[0]]<MAX_COUNT)
								patternCounts[-KWTree[currentNodeID].children[0]]++;
						}
						currentPositionInInput++;
						found=1;
					}
					else
					{
						suffixLinkID= KWTree[currentNodeID].suffixLinkID;  //follow up
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

	return EXIT_SUCCESS;
}


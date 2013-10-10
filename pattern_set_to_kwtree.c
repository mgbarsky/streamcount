#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "KeywordTree.h"

//prepriocesses k-mers created from one patterns file inrto a keyword tree
int preprocessPatternSet(KWTreeBuildingManager *manager, char *patternsFileName, int availableRamMB, int64_t *totalUniquePatterns)
{
	//find out what is the maximum size of keyword tree we can build in memory 
	//- so limiting the total length of pattern set
	char currentFileName[MAX_PATH_LENGTH];
	char ** patterns;
	int i;
	int64_t totalPatterns=0LL;
	uint64_t memBytes = availableRamMB *1000000;	
	FILE *patternsInputFP;
	FILE *outputFP;

	uint64_t maxSetSize = memBytes/sizeof(KWTNode);
	manager->maxSetSize = MIN(maxSetSize,MAX_SET_SIZE);

	//init fields
	manager->treeSlotsNum =1; //next available slot, the root is in a slot 0	
	
	//here we want to know how much memory to allocate for a kw tree
	//open patterns file for reading
	if(!(patternsInputFP= fopen ( patternsFileName , "r" )))
	{
		printf("Could not open file \"%s\" for reading patterns \n", patternsFileName);
		return 1;
	}	
		
	if(collectPatternsStats(patternsInputFP,manager->k,&totalPatterns)!=0)
	{
		fclose(patternsInputFP);
		return 1;
	}
	
	printf("Total %ld k-mers of size %d each\n", totalPatterns,(manager->k)); 
	
	//check we have sufficient memory to hold KWtree
	if(totalPatterns*(manager->k+1) < manager->maxSetSize)
		manager->totalSetSize = totalPatterns*(manager->k+1);
	else
	{
		printf("Pattern set can not be preprocessed in the amount of the available memory:\n");
		printf("We can hold maximum %d keyword tree nodes, but the tree may require %ld nodes\n",
				manager->maxSetSize,totalPatterns*(manager->k+1));
		return 1;
	}
	
	//allocate an array to hold all patterns
	if(!(patterns =(char **)calloc(totalPatterns+1,sizeof (char *)))) //=1 because we will store unique patterns starting from position 1, and what if all totalPatterns are unique
	{
		printf("Unable to allocate memory for patterns array.\n");
		return 1;
	}
	for(i=0;i<totalPatterns+1;i++)
	{
		if(!(patterns[i] =(char *)calloc(manager->k+1,sizeof (char *))))
		{
			printf("Unable to allocate memory for pattern %d.\n",i);
			return 1;
		}
	}
	
	rewind(patternsInputFP);
	//fill-in array of patterns
	if(fillPatternsArray(patternsInputFP, patterns, &totalPatterns,manager->k)!=0)
		return 1;
	//close input file
	fclose(patternsInputFP);	
	
	//allocate memory for KWtree
	if(!(manager->KWtree =(KWTNode *)calloc(manager->totalSetSize,sizeof (KWTNode))))
	{
		printf("Unable to allocate memory for keyword tree.\n");
		return 1;
	}
	
	//pre-process patterns into a keyword tree -there duplicate k-mers get removed	
	if(buildKeywordTree(manager, patterns,&totalPatterns, totalUniquePatterns)!=0)
		return 1;

	printf("BuildKeywordTree complete\n");

	//write unique sorted patterns to file - mapping from line number to actual pattern 
	//open file for writing
	sprintf(currentFileName, "%s_%d-mers_MAPPINGS",  patternsFileName,manager->k);
	if(!(outputFP= fopen ( currentFileName , "w" )))
	{
		printf("Could not open file \"%s\" for writing unique patterns \n", currentFileName);
		return 1;
	}
	
	//write patterns to file
	for(i=0;i<*totalUniquePatterns;i++)
	{
		fprintf(outputFP,"%s\n",patterns[i]);
	}
	fclose(	outputFP);
	
	//free patterns array - dont need it anymore 
	for(i=0;i<totalPatterns;i++)
	{
		free(patterns[i]); 
	}
	free(patterns);	

	return 0;
}


//Reads though input patterns file and collects the estimated number of k-mers - into *totalPatterns - to allocate memory
int collectPatternsStats(FILE *inputFP, int k, int64_t *totalPatterns)
{
	char currentLine[MAX_CHARS_PER_LINE];
	*totalPatterns=0;
	while( fgets (currentLine, MAX_CHARS_PER_LINE-10, inputFP)!=NULL ) 
	{
		int lineLen = strlen(currentLine);
		*totalPatterns+=(lenValidChars(currentLine,lineLen)-k+1);
	}	
	return 0;
}

//reads input patterns file again and creates non-unique k-mers from each input line
int fillPatternsArray(FILE *inputFP, char **patterns, int64_t *totalPatterns, int k)
{
	char currentLine[MAX_CHARS_PER_LINE];
	int i;
	int patternsCounter=0;
	
	while( fgets (currentLine, MAX_CHARS_PER_LINE-10, inputFP)!=NULL ) 
	{
		int lineLen = lenValidChars(currentLine,strlen(currentLine));
		
		for(i=0;i<lineLen-k+1;i++)
		{
			memcpy(patterns[patternsCounter],&currentLine[i],k);
			patterns[patternsCounter][k]='\0';
			patternsCounter++;
		}			
	}	
	
	//sanity check
	if(patternsCounter!=*totalPatterns)
	{
		printf("Unexpected error: estimated number of patterns is %ld but added %d\n",*totalPatterns,patternsCounter);
		return 1;
	}
	return 0;
}


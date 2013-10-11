#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "KeywordTree.h"

//Builds keyword tree with suffix links for each k-mer
//obtained from a set of patterns extracted from file globalArgs.patternFileName
//Converting each line of the input file (line length cannot exceed MAX_CHARS_PER_LINE) into k-mer substrings
int 
preprocessPatternsIntoKeywordTree(GlobalArgs *globalArgs)
{
	char patternFileName [MAX_PATH_LENGTH]; 
	KWTreeBuildingManager manager;
	int written;
	char currentFileName[MAX_PATH_LENGTH];
	FILE *outputFP;
	KWTreeInfo info[1];
	int64_t totalUniquePatterns=0;
	
	int memoryMB;
	
	
	sprintf(patternFileName,"%s", globalArgs->patternFileName);	
	manager.k = globalArgs->k;	
	memoryMB = globalArgs->memoryInMB;
	
	printf("Processing patterns file %s into a set of unique %d-mers using %d MB of RAM\n",patternFileName,manager.k,memoryMB);

	if(preprocessPatternSet(&manager, patternFileName, memoryMB, &totalUniquePatterns)!=0)
	{
		printf("Error building KWtree\n");
		return 1;
	}	
	printf("Preprocessing into KWtree complete\n");

	//now manager contains an important information about pattern set and the kwtree itself
	//we serialize this information to store for future processing with each file
	info[0].treeSlotsNum = manager.treeSlotsNum;
	info[0].k=manager.k;
	info[0].totalPatterns = totalUniquePatterns;

	sprintf(currentFileName, "%s_%d-mers_KWTREE_INFO", patternFileName, manager.k);

	if(!(outputFP= fopen ( currentFileName , "wb" )))
	{
		printf("Could not open file \"%s\" for writing KWtree info \n", currentFileName);
		return 1;
	}

	if(fwrite (&info[0] , sizeof(KWTreeInfo),1, outputFP)!=1)
	{
		printf("Could not save KWtree info \n");
		return 1;
	}

	printf("Saved KWtree info to file %s: totalNodes = %d, k=%d, totalPatterns=%d\n",
			currentFileName,info[0].treeSlotsNum,info[0].k,info[0].totalPatterns);
	fclose(outputFP);
	
	//serialize KWtree to disk
	sprintf(currentFileName, "%s_%d-mers_KWTREE", patternFileName, manager.k);
	if(!(outputFP= fopen ( currentFileName , "wb" )))
	{
		printf("Could not open file \"%s\" for writing KWtree NODES \n", currentFileName);
		return 1;
	}

	written = fwrite (&(manager.KWtree[0]) , sizeof(KWTNode), manager.treeSlotsNum, outputFP);

	if(written!=manager.treeSlotsNum)
	{
		printf("Not all KWtree nodes have been written: wanted to write %d, and wrote %d \n",manager.treeSlotsNum,written);
		return 1;
	}
	fclose(outputFP);	

	

	//free memory
	freeMemoryAfterKWtBuild (&manager);

	printf("Happy End\n");
	return 0;
}

int freeMemoryAfterKWtBuild (KWTreeBuildingManager* manager)
{
	if(manager->KWtree)
		free(manager->KWtree);
	return 0;
}

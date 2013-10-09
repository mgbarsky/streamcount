#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "KeywordTree.h"

//Builds keyword tree with suffix links for each k-mer
//obtained from a set of patterns specified in <patterns file name fullpath>
//Converting each line of the input file (line length cannot exceed MAX_CHARS_PER_LINE) into k-mer substrings
//Parameters: full path and name of the patterns file, value of K, available memory in MB
int 
main(int argc, char *argv[])
{
	char patternFileName [MAX_PATH_LENGTH]; 
	KWTreeBuildingManager manager;
	int written;
	char currentFileName[MAX_PATH_LENGTH];
	FILE *outputFP;
	KWTreeInfo info[1];
	int64_t totalUniquePatterns=0;
	KWTNode *tree;
	int memoryMB;

	if(argc<4)
	{
		printf("./buildkwtree <patterns file name fullpath> <k-mer size> <memory_megabytes_available_for_kwtree> \n");
		return 1;
	}
	
	sprintf(patternFileName,"%s", argv[1]);	
	manager.k = atoi(argv[2]);	
	memoryMB = atoi(argv[3]);
	
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

	//here we are going to test if the KWTree is the same when deserialized
	if(DEBUG_SERIALIZATION)
	{
		//read info
		KWTreeInfo newinfo[1];
		FILE *inputFP;
		int read;
		//read info from file
		sprintf(currentFileName, "%s_%d-mers_KWTREE_INFO", patternFileName, manager.k);
		if(!(inputFP= fopen ( currentFileName , "rb" )))
		{
			printf("Could not open file \"%s\" for reading saved KWtree info \n", currentFileName);
			return 1;
		}
		
		//read info from file
		read= fread (newinfo,sizeof(KWTreeInfo),1,inputFP);

		if(read!=1)
		{
			printf("Error reading KWTinfo from file\n");
			return 1;
		}
		fclose(inputFP);
		//print new tree parameters read from file
		printf("Read from file %s the following parameters: : totalNodes = %d, k=%d, totalPatterns=%d\n",
			currentFileName,newinfo[0].treeSlotsNum,newinfo[0].k,newinfo[0].totalPatterns);
		

		//allocate memory for KWtree		
		if(!(tree =(KWTNode *)calloc(newinfo[0].treeSlotsNum,sizeof (KWTNode))))
		{
			printf("Unable to allocate memory to reload keyword tree.\n");
			return 1;
		}
		
		//open file for reading a tree into RAM
		sprintf(currentFileName, "%s_%d-mers_KWTREE", patternFileName, newinfo[0].k);
		if(!(inputFP= fopen ( currentFileName , "rb" )))
		{
			printf("Could not open file \"%s\" for reading saved KWtree \n", currentFileName);
			return 1;
		}
		
		//read tree from file
		read= fread (tree,sizeof(KWTNode),newinfo[0].treeSlotsNum,inputFP);
		if(read!=newinfo[0].treeSlotsNum)
		{
			printf("Error reading KWT tree from file: wanted to read %d nodes but read %d\n", newinfo[0].treeSlotsNum, read);
			return 1;
		}
		free(tree);
	}

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

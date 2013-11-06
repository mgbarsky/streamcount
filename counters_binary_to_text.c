#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include "KeywordTree.h"


int 
main(int argc, char *argv[])
{
	char patternFileName [MAX_PATH_LENGTH];
	int i,w, k;
	
	char currentFileName [MAX_PATH_LENGTH];

	KWTreeInfo info[1];
	INT totalCounters;
	FILE *inputFP, *outputFP;
	int read;
	INT *counters;

	if(argc<3)
	{
		printf("./countstotext <patterns file name fullpath> <size of k-mers> <input file name(s) full path>\n");
		return 1;
	}
	
	snprintf(patternFileName, MAX_PATH_LENGTH, "%s", argv[1]);
	
	k=atoi(argv[2]);

	//try to read KWtree info file to find out the size of the counters array
	snprintf(currentFileName, MAX_PATH_LENGTH, "%s_%d-mers_KWTREE_INFO", patternFileName, k);
	if(!(inputFP= fopen ( currentFileName , "rb" )))
	{
		printf("Could not open file \"%s\" for reading saved KWtree info \n", currentFileName);
		return 1;
	}
		
	read= fread (info,sizeof(KWTreeInfo),1,inputFP);
	if(read!=1)
	{
		fclose(inputFP);
		printf("Error reading KWTinfo from file\n");
		return 1;
	}
	fclose(inputFP);

	//print how many leaves in the tree
	totalCounters=info[0].totalLeaves;
	printf("Total counters=%d\n",totalCounters);

	//allocate an array of counters
	if(!(counters =(INT *)calloc(totalCounters,sizeof (INT))))
	{
		printf("Unable to allocate memory for counters array.\n");
		return 1;
	}

	//go in the loop through remaining arguments which are the names of the processed files
	for(i=3; i<argc; i++)
	{
		snprintf(currentFileName, MAX_PATH_LENGTH, "%s_%d-mers_COUNTS", argv[i], k);
		if(!(inputFP= fopen ( currentFileName , "rb" )))
		{
			printf("Could not open file \"%s\" for reading binary counters \n", currentFileName);
			return 1;
		}
		
		read=fread (&(counters[0]),sizeof(INT),totalCounters,inputFP);
		if(read!=totalCounters)
		{
			printf("Error reading counts back from file: wanted to read %d but read %d\n", totalCounters, read);
			return 1;
		}
		fclose(inputFP);

		//open text file for writing
		snprintf(currentFileName, MAX_PATH_LENGTH, "%s_%d-mers_COUNTS.txt", argv[i], k);
		if(!(outputFP= fopen ( currentFileName , "w" )))
		{
			printf("Could not open file \"%s\" for whiting textual counters \n", currentFileName);
			return 1;
		}

		for(w=0;w<totalCounters;w++)
		{
			INT count = counters[w];
			fprintf(outputFP,"%d\n", count);
		}
		fclose(outputFP);
	}
	return 0;
}


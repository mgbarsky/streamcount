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
	int  i,k;
	
	char currentFileName [MAX_PATH_LENGTH];

	
	INT totalPatterns;
	int bytes;
	FILE *inputFP, *outputFP;
	int read;
	KmerInfo *patternsInfo;
	
	if(argc<3)
	{
		printf("./patternsstotext <patterns file name fullpath> <size of k-mers> \n");
		return 1;
	}
	
	snprintf(patternFileName, MAX_PATH_LENGTH, "%s", argv[1]);
	
	k=atoi(argv[2]);

	//try to read k-mers info 
	snprintf(currentFileName, MAX_PATH_LENGTH, "%s_%d-mers_MAPPINGS",  patternFileName, k);
	if(!(inputFP= fopen ( currentFileName , "rb" )))
	{
		printf("Could not open file \"%s\" for reading k-mer mapping \n", currentFileName);
		return 1;
	}
	
	fseek (inputFP, 0, SEEK_END);   // non-portable
	bytes = ftell (inputFP);
	printf("%d bytes\n",bytes);
    	totalPatterns=bytes/sizeof(KmerInfo);
	
	printf("Total k-mers=%d\n",totalPatterns);

	//allocate an array of patterns info
	if(!(patternsInfo =(KmerInfo *)calloc(totalPatterns,sizeof (KmerInfo))))
	{
		printf("Unable to allocate memory for patterns info array.\n");
		return 1;
	}
	
	rewind(inputFP);

	read= fread (patternsInfo,sizeof(KmerInfo),totalPatterns,inputFP);
	if(read!=totalPatterns)
	{
		fclose(inputFP);
		printf("Error reading pattern info from file\n");
		return 1;
	}
	fclose(inputFP);
	
	//open text file for writing
	snprintf(currentFileName, MAX_PATH_LENGTH, "%s_%d-mers_MAPPINGS.txt",  patternFileName, k);
	if(!(outputFP= fopen ( currentFileName , "w" )))
	{
		printf("Could not open file \"%s\" for whriting textual k-mer info \n", currentFileName);
		return 1;
	}
	
	//go in a loop through k-mers and print their info to file
	for(i=0; i<totalPatterns; i++)
	{
		fprintf(outputFP,"%d\t%d\t%d\t%d\t%d\n", patternsInfo[i].lineNumber, patternsInfo[i].startPosInLine, patternsInfo[i].counterID,patternsInfo[i].rcCounterID, patternsInfo[i].repeated);		
	}
	fclose(outputFP);
	return 0;
}


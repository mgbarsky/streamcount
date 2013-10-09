#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include "KeywordTree.h"


int streamAndCountOneFile(KWTCounterManager *manager)
{
	FILE *inputFP;
	char inputFileName[MAX_PATH_LENGTH];
	char currentLine[MAX_CHARS_PER_LINE];
	INT currentLineAsINT[MAX_CHARS_PER_LINE];

	sprintf(inputFileName, "%s%d", manager->inputFileNamePrefix, manager->currentFileID);

	//open file to read lines
	if(!(inputFP= fopen ( inputFileName , "r" )))
	{
		printf("Could not open input file \"%s\" for reading\n", inputFileName);
		return 1;
	}

	while( fgets (currentLine, MAX_CHARS_PER_LINE-10, inputFP)!=NULL ) 
	{
		int lineLen = validCharsToIntArray(currentLine,strlen(currentLine),&currentLineAsINT[0]);
		
		if(streamOneString(manager,currentLineAsINT,lineLen)!=0)
		{
			fclose(inputFP);
			return 1;
		}	
	}	
	fclose(inputFP);
	return 0;	
}

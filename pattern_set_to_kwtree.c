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
	KmerInfo *patternsInfo;
	int i;
	int written;
	int64_t totalPatterns=0LL;
	uint64_t memBytes = availableRamMB *1000000;	
	FILE *patternsInputFP;
	FILE *outputFP;

	int totalKMers = 0;

	uint64_t maxSetSize = memBytes/sizeof(KWTNode);
	manager->maxSetSize = MIN(maxSetSize,MAX_SET_SIZE);
	
	

	//init fields of a tree
	manager->treeSlotsNum =1; //next available slot, the root is in a slot 0	
	manager->treeLeavesNum=1; //we start from -1

	//here we want to know how much memory to allocate for a kw tree
	//open patterns file for reading
	if(!(patternsInputFP= fopen ( patternsFileName , "r" )))
	{
		printf("Could not open file \"%s\" for reading patterns \n", patternsFileName);
		return 1;
	}	
		
	if(collectPatternsStats(patternsInputFP,manager->k,manager->inputType, &totalPatterns)!=0)
	{
		fclose(patternsInputFP);
		return 1;
	}
	
	manager->totalPatterns = totalPatterns;
	printf("Estimated total %ld non-rc k-mers of size %d each\n", totalPatterns,(manager->k)); 
	
	if(totalPatterns == 0)
	{
		printf("No valid characters found in the pattern input file, or the k-mer size exceeds the line length\n");
		return 1;
	}
	//check we have sufficient memory to hold KWtree

	totalKMers = (manager->includeReverseComplement)?2*totalPatterns:totalPatterns;

	if(totalKMers*(manager->k+1) < manager->maxSetSize)
		manager->totalSetSize = totalKMers*(manager->k+1);
	else
	{
		printf("Pattern set can not be preprocessed in the amount of the available memory:\n");
		printf("We can hold maximum %d keyword tree nodes, but the tree may require %d nodes\n",
				manager->maxSetSize,totalKMers*(manager->k+1));
		return 1;
	}
	
	//allocate an array to hold all patterns and their info
	if(!(patterns =(char **)calloc(totalPatterns,sizeof (char *)))) 
	{
		printf("Unable to allocate memory for patterns array.\n");
		return 1;
	}

	for(i=0;i<totalPatterns;i++)
	{
		if(!(patterns[i] =(char *)calloc(manager->k+1,sizeof (char *))))
		{
			printf("Unable to allocate memory for pattern %d.\n",i);
			return 1;
		}
	}
	
	if(!(patternsInfo =(KmerInfo *)calloc(totalPatterns,sizeof (KmerInfo)))) 
	{
		printf("Unable to allocate memory for patterns info array.\n");
		return 1;
	}

	rewind(patternsInputFP);
	//fill-in array of patterns
	if(fillPatternsArray(patternsInputFP, patterns, patternsInfo,  &totalPatterns, manager->k, manager->inputType)!=0)
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
	if(buildKeywordTree(manager, patterns, patternsInfo, &totalPatterns, totalUniquePatterns)!=0)
		return 1;

	printf("BuildKeywordTree complete. Contains %ld leaves\n", *totalUniquePatterns);

	//write patterns info to file - mapping from line number and position to actual leaf - position in the array of counters 
	//open file for writing
	snprintf(currentFileName, MAX_PATH_LENGTH, "%s_%d-mers_MAPPINGS",  patternsFileName,manager->k);
	if(!(outputFP= fopen ( currentFileName , "wb" )))
	{
		printf("Could not open file \"%s\" for writing k-mer mapping \n", currentFileName);
		return 1;
	}
	
	//write patterns info to file
	printf("Writing %ld patterns info to file %s\n",totalPatterns,currentFileName); 
	written = fwrite (&(patternsInfo[0]) , sizeof(KmerInfo), totalPatterns, outputFP);

	if(written!=totalPatterns)
	{
		printf("Not all pattern info was written: wanted to write %ld, and wrote %d \n",totalPatterns,written);
		return 1;
	}
	fclose(	outputFP);
	
	if(DEBUG_KMERS)
	{
		printf("******Summary of k-mers collected from the input*******\n");
		for(i=0;i<totalPatterns;i++)
		{
			printf("%d. pattern %s, counterID=%d, rcCounterID=%d, line=%d, pos=%d, repeated=%d\n",i,patterns[i],
								patternsInfo[i].counterID,patternsInfo[i].rcCounterID,
								patternsInfo[i].lineNumber,patternsInfo[i].startPosInLine,patternsInfo[i].repeated);

		}
	}
	//free patterns array - dont need it anymore 
	for(i=0;i<totalPatterns;i++)
	{
		free(patterns[i]); 
	}
	free(patterns);	
	free(patternsInfo);

	return 0;
}


//Reads though input patterns file and collects the estimated number of k-mers - into *totalPatterns - to allocate memory
int collectPatternsStats(FILE *inputFP, int k, int inputType, int64_t *totalPatterns)
{
	char currentLine[MAX_CHARS_PER_LINE];
	*totalPatterns=0;
	
	while( fgets (currentLine, MAX_CHARS_PER_LINE-10, inputFP)!=NULL ) 
	{
		int lineLen = strlen(currentLine);
		switch( inputType ) 
		{
            		case INPUT_LINES:  //lines
                		*totalPatterns+=MAX(0,(lenValidChars(currentLine,lineLen)-k+1));
                		break;

			case INPUT_FILE:  //entire file
				*totalPatterns+=MAX(0,(lenValidChars(currentLine,lineLen)));
				break;
			/*case INPUT_SNIPS: //snips
				*totalPatterns+=numberKmersFromSnips(currentLine,lineLen,k);
				break;*/
			default:
				printf("UNEXPECTED INPUT TYPE %d\n",inputType);
				return 1;
		}		
	}	
	return 0;
}

//reads input patterns file again and creates non-unique k-mers from each input line
int fillPatternsArray(FILE *inputFP, char **patterns, KmerInfo *patternsInfo, int64_t *totalPatterns, int k, int inputType)
{
	char currentLine[MAX_CHARS_PER_LINE];
	char previousLine [MAX_CHARS_PER_LINE];
	int i,j;
	int patternsCounter=0;
	int lineCounter = 0;
	int prevLineLen =0;

	while( fgets (currentLine, MAX_CHARS_PER_LINE-10, inputFP)!=NULL ) 
	{
		if(inputType == INPUT_LINES)
		{
			int lineLen = lenValidChars(currentLine,strlen(currentLine));
			if(lineLen>=k)
			{
				for(i=0;i<lineLen-k+1;i++)
				{
					memcpy(patterns[patternsCounter],&currentLine[i],k);
					patterns[patternsCounter][k]='\0';
					patternsInfo[patternsCounter].lineNumber = lineCounter;
					patternsInfo[patternsCounter].startPosInLine = i;
					patternsInfo[patternsCounter].repeated = 0; //initialize repetitions count
					patternsCounter++;

					if(patternsCounter > *totalPatterns)
					{
						printf("Unexpected error while extracting k-mers: too many k-mers found\n");
						return 1;
					}
				}
				
			}
			lineCounter++;			
		}
		else if (inputType == INPUT_FILE)
		{
			int lineLen = lenValidChars(currentLine,strlen(currentLine));
			
			if(lineLen>3)
			{
				//extract k-mers which start at k-1 positions of a previous line
				if(prevLineLen>0) //there is a previous line - we need to concatenate suffix of this line with prefix of the current line
				{
					for(i=prevLineLen - (k-1); i< prevLineLen; i++)
					{
						int prefixLen = prevLineLen -i ;
						//check that there is enough suffix len in the current line to concatenate
						if(lineLen >= k-prefixLen)
						{
							//copy prefix
							memcpy(patterns[patternsCounter],&previousLine[i],prefixLen);

							//copy suffix
							for(j=prefixLen; j < k; j++)
								patterns[patternsCounter][j] = currentLine[j - prefixLen];

							patterns[patternsCounter][k]='\0';
							patternsInfo[patternsCounter].lineNumber = lineCounter-1;  //starts in the previous line
							patternsInfo[patternsCounter].startPosInLine = i;
							patternsInfo[patternsCounter].repeated = 0; //initialize repetitions count
							patternsCounter++;

							if(patternsCounter > *totalPatterns)
							{
								printf("Unexpected error while extracting k-mers: too many k-mers found\n");
								return 1;
							}
						}
					}
				}
			
				//process current line, if it has enough length, upto k-1 last characters
				if(lineLen>=k)
				{
					for(i=0;i<lineLen-k+1;i++)
					{
						memcpy(patterns[patternsCounter],&currentLine[i],k);
						patterns[patternsCounter][k]='\0';
						patternsInfo[patternsCounter].lineNumber = lineCounter;
						patternsInfo[patternsCounter].startPosInLine = i;
						patternsInfo[patternsCounter].repeated = 0; //initialize repetitions count
						patternsCounter++;

						if(patternsCounter > *totalPatterns)
						{
							printf("Unexpected error while extracting k-mers: too many k-mers found\n");
							return 1;
						}
					}
					
					//keep this line to be concatenated with the next line
					memcpy(&previousLine[0],&currentLine[0],lineLen);
					previousLine[lineLen]='\0';
					prevLineLen=lineLen;
				}				
			}
			else  //either empty line or description line
			{
				prevLineLen=0;
			}
			lineCounter++;
		}

		/*else if(inputType == INPUT_SNIPS)
		{
			int lineLen = lenValidChars(currentLine,strlen(currentLine));
			
			//first, generate all k-mers from line without snips
			for(i=0;i<lineLen-k+1;i++)
			{
				memcpy(patterns[patternsCounter],&currentLine[i],k);
				patterns[patternsCounter][k]='\0';
				patternsInfo[patternsCounter].lineNumber = lineCounter;
				patternsInfo[patternsCounter].startPosInLine = i;
				patternsInfo[patternsCounter].repeated = 0; //initialize repetitions count
				patternsCounter++;

				if(patternsCounter > *totalPatterns)
				{
					printf("Unexpected error while extracting k-mers: too many k-mers found\n");
					return 1;
				}
			}

			//now we need to add all k-mers (k-1) characters to the left of a specified polymorphic position, if we replace one nucleotyde by its polymorphic base
			i=lineLen;
			
			int middle=0;
			int snipLen;
			int snipVal;
			while(currentLine[i]==',' && validSnip(&(currentLine[i+1]),lineLen,&middle, &snipLen, &snipVal))
			{
				for(j=MAX(0,middle - (k-1)); j <= MIN(middle,lineLen - k); j++)
				{
					memcpy(patterns[patternsCounter],&currentLine[j],k);
					patterns[patternsCounter][middle-j] = (char)snipVal;
					patterns[patternsCounter][k]='\0';
					patternsInfo[patternsCounter].lineNumber = lineCounter;
					patternsInfo[patternsCounter].startPosInLine = -(i+1);
					patternsInfo[patternsCounter].repeated = 0; //initialize repetitions count
					patternsCounter++;

					if(patternsCounter > *totalPatterns)
					{
						printf("Unexpected error while extracting k-mers: too many k-mers found\n");
						return 1;
					}
				}
				i+=snipLen;
			}
		}*/			
	}
	
	*totalPatterns = patternsCounter;
	return 0;
}


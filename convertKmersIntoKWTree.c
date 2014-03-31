/*
This will extract all possible k-mers from the input file
and convert them into a disk-resident keyword tree
saving information about k-mers - their mapping to input lines
*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include "General.h"
#include "StreamCount.h"


int convertAllKmersIntoKWTreeReturnTree (FILE *kmersFP, int inputType, INT k, int includeRC, INT memoryMB, KWTreeBuildingManager *manager)
{
	INT estimatedNumberOfKmers=0;
	INT estimatedTotalNumberOfKmers; //twice estimatedNumberOfKmers if include rc	
	INT i;
	
    fprintf(stderr,"________________________________________\n");
	fprintf(stderr,"Started processing input %ld-mers into an index using %ld MB of RAM\n",
		(long)manager->k,(long)memoryMB);
    fprintf(stderr,"****************************\n");

	//this will read input file and estimate how many k-mers it is possible to extract - return into estimatedNumberOfKmers
	if(collectKmerInputStats(kmersFP,manager->k,manager->inputType, &estimatedNumberOfKmers)!=EXIT_SUCCESS)
	{
		fclose(kmersFP);
		return EXIT_FAILURE;
	}

	if(DEBUG_KMERS_EXTRACTION) fprintf(stderr,"Estimated total %ld non-rc k-mers of size %ld each\n", (long)estimatedNumberOfKmers,(long)(manager->k)); 
	
	if(estimatedNumberOfKmers == 0)
	{
		fprintf(stderr,"No valid characters found in the pattern input file, or the k-mer size exceeds the line length\n");
		return EXIT_FAILURE;
	}	
	
	//check if we have sufficient memory to hold KWtree
	estimatedTotalNumberOfKmers = (manager->includeReverseComplement)?2*estimatedNumberOfKmers:estimatedNumberOfKmers;
	manager->estimatedNumberOfLeaves = estimatedTotalNumberOfKmers+1;

	if(estimatedTotalNumberOfKmers*(manager->k+1) < manager->maxSetSize && estimatedTotalNumberOfKmers*(manager->k+1) < MAX_INT)
		manager->estimatedNumberOfKWTreeNodes = estimatedTotalNumberOfKmers*(manager->k+1);
	else
	{
		fprintf(stderr,"Pattern set can not be preprocessed in the amount of the available memory:\n");
		fprintf(stderr,"We can hold maximum %ld keyword tree nodes, but the tree may require %ld nodes\n",
				(long)manager->maxSetSize,(long)estimatedTotalNumberOfKmers*(manager->k+1));
		return EXIT_FAILURE;
	}
	
	//allocate an array to hold all kmers and their info
	if(!(manager->kmers =(char **)calloc(estimatedNumberOfKmers,sizeof (char *)))) 
	{
		fprintf(stderr,"Unable to allocate memory for kmers array.\n");
		return EXIT_FAILURE;
	}

	for(i=0;i<estimatedNumberOfKmers;i++)
	{
		if(!(manager->kmers[i] =(char *)calloc(manager->k+1,sizeof (char *))))
		{
			fprintf(stderr,"Unable to allocate memory for kmer %ld.\n",(long)i);
			return EXIT_FAILURE;
		}
	}
	
	if(!(manager->kmersInfo =(KmerInfo *)calloc(estimatedNumberOfKmers,sizeof (KmerInfo)))) 
	{
		fprintf(stderr,"Unable to allocate memory for kmers info array.\n");
		return EXIT_FAILURE;
	}

    //remember estimatedNumberOfKmers to free these arrays later
    manager->estimatedNumberOfKmers = estimatedNumberOfKmers;

	rewind(kmersFP); 

	//after this, the actual number of valid k-mers to process is in manager->originalNumberOfKmers
	if(fillKmersArrayAndInfo(manager, kmersFP, manager->kmers, manager->kmersInfo, estimatedNumberOfKmers)!=EXIT_SUCCESS)
	{
		fprintf(stderr,"Error generating kmers array\n");
		fclose(kmersFP);
		return EXIT_FAILURE;
	}
	
	if (DEBUG_KMERS_EXTRACTION) fprintf(stderr,"Preprocessing of the input file into an array of kmers complete\n");

	//allocate memory for KWtree
	INT maxNumberOfNodes = (manager->maxNumberOfLeaves)*(manager->k)+1; //each of k chars can be potentially a node plus a root node
	if(!(manager->KWtree =(KWTNode *)calloc(maxNumberOfNodes,sizeof (KWTNode))))
	{
		fprintf(stderr,"Unable to allocate memory for keyword tree.\n");
		return EXIT_FAILURE;
	}
    
	if (DEBUG_KWTREE)
    {
        INT allMemInBytes=maxNumberOfNodes*sizeof (KWTNode);
        fprintf(stderr,"Allocated memory %ld to hold %ld KW tree nodes with max %ld number of leaves \n",(long)allMemInBytes,(long) maxNumberOfNodes,(long)manager->maxNumberOfLeaves);
    }

	//pre-process patterns into a keyword tree -there duplicate k-mers get removed
	if(buildKeywordTree(manager, manager->kmers, manager->kmersInfo)!=EXIT_SUCCESS)
		return EXIT_FAILURE;

	if(DEBUG_KWTREE) fprintf(stderr,"BuildKeywordTree complete. totalNodes = %ld, k=%ld, totalLeaves=%ld\n", (long)manager->actualNumberOfKWTreeNodes,(long)manager->k,(long)manager->actualNumberOfLeaves);	
	
    fprintf(stderr,"***********************************\n");
	fprintf(stderr,"Happy end for KWTree\n");
    fprintf(stderr,"________________________________________\n\n");
	return EXIT_SUCCESS;
}


int collectKmerInputStats(FILE *inputFP, INT k, int inputType, INT *estimatedNumberOfKmers)
{
	char currentLine[MAX_CHARS_PER_LINE];
	INT kmersCount=0;
	
	while( fgets (currentLine, MAX_CHARS_PER_LINE-10, inputFP)!=NULL ) 
	{
		int lineLen = strlen(currentLine);
		switch( inputType ) 
		{
            case KMERS_FROM_LINES:  //lines
                kmersCount+=MAX(0,(lenValidChars(currentLine,lineLen)-k+1));
                break;

			case KMERS_FROM_FILE:  //entire file
				kmersCount+=MAX(0,(lenValidChars(currentLine,lineLen)));
				break;
			
			default:
				fprintf(stderr,"UNEXPECTED INPUT TYPE %d\n", inputType);
				return EXIT_FAILURE;
		}		
	}	
	*estimatedNumberOfKmers = kmersCount;
	return EXIT_SUCCESS;
}

int fillKmersArrayAndInfo(KWTreeBuildingManager* manager, FILE *inputFP, 
	char **kmers, KmerInfo *kmersInfo, INT maxPossibleNumberOfKmers)
{
	char currentLine[MAX_CHARS_PER_LINE];
	char previousLine [MAX_CHARS_PER_LINE];
	int i,j;
	int kmersCounter=0;
	int lineCounter = 0;
	int prevLineLen =0;

    if (DEBUG_KWTREE) fprintf(stderr,"maxPossibleNumberOfKmers=%ld\n",(long)maxPossibleNumberOfKmers);
	while( fgets (currentLine, MAX_CHARS_PER_LINE-10, inputFP)!=NULL ) 
	{
		if(manager->inputType == KMERS_FROM_LINES)
		{
			int lineLen = lenValidChars(currentLine,strlen(currentLine));
			if(lineLen>=manager->k)
			{
				for(i=0;i<lineLen-manager->k+1;i++)
				{
					memcpy(kmers[kmersCounter],&currentLine[i],manager->k);
					kmers[kmersCounter][manager->k]='\0';
					kmersInfo[kmersCounter].lineNumber = lineCounter;
					kmersInfo[kmersCounter].startPosInLine = i;
					kmersInfo[kmersCounter].repeated = 0; //initialize repetitions count
					kmersCounter++;

					if(kmersCounter > maxPossibleNumberOfKmers) //to avoid memory overflow
					{
						fprintf(stderr,"Unexpected error while extracting k-mers: too many k-mers found\n");
						return EXIT_FAILURE;
					}
				}
				
			}
            else
            {
                fprintf(stderr,"Line %d contains non-DNA character or is less than k=%ld: %s\n",lineCounter,(long)manager->k,currentLine);
                return EXIT_FAILURE;
            }
			lineCounter++;			
		}
		else if (manager->inputType == KMERS_FROM_FILE)
		{
			int lineLen = lenValidChars(currentLine,strlen(currentLine));
			
			if(lineLen>3)
			{
				//extract k-mers which start at k-1 positions of a previous line
				if(prevLineLen>0) //there is a previous line - we need to concatenate suffix of this line with prefix of the current line
				{
					for(i=prevLineLen - (manager->k-1); i< prevLineLen; i++)
					{
						int prefixLen = prevLineLen -i ;
						//check that there is enough suffix len in the current line to concatenate
						if(lineLen >= manager->k-prefixLen)
						{
							//copy prefix
							memcpy(kmers[kmersCounter],&previousLine[i],prefixLen);

							//copy suffix
							for(j=prefixLen; j < manager->k; j++)
								kmers[kmersCounter][j] = currentLine[j - prefixLen];

							kmers[kmersCounter][manager->k]='\0';
							kmersInfo[kmersCounter].lineNumber = lineCounter-1;  //starts in the previous line
							kmersInfo[kmersCounter].startPosInLine = i;
							kmersInfo[kmersCounter].repeated = 0; //initialize repetitions count
							kmersCounter++;

							if(kmersCounter > maxPossibleNumberOfKmers) //to avoid memory overflow
							{
								fprintf(stderr,"Unexpected error while extracting k-mers: too many k-mers found\n");
								return EXIT_FAILURE;
							}
						}
					}
				}
			
				//process current line, if it has enough length, upto k-1 last characters
				if(lineLen>=manager->k)
				{
					for(i=0;i<lineLen-manager->k+1;i++)
					{
						memcpy(kmers[kmersCounter],&currentLine[i],manager->k);
						kmers[kmersCounter][manager->k]='\0';
						kmersInfo[kmersCounter].lineNumber = lineCounter;
						kmersInfo[kmersCounter].startPosInLine = i;
						kmersInfo[kmersCounter].repeated = 0; //initialize repetitions count
						kmersCounter++;

						if(kmersCounter > maxPossibleNumberOfKmers)
						{
							fprintf(stderr,"Unexpected error while extracting k-mers: too many k-mers found\n");
							return EXIT_FAILURE;
						}
					}
					
					//keep this line to be concatenated with the next line
					memcpy(&previousLine[0],&currentLine[0],lineLen);
					previousLine[lineLen]='\0';
					prevLineLen=lineLen;
				}
                else
                {
                    fprintf(stderr,"Line %d contains non-DNA character or is less than k=%ld: %s\n",lineCounter,(long)manager->k,currentLine);
                    return EXIT_FAILURE;
                }				
			}
			else  //either empty line or description line
			{
				prevLineLen=0;
			}
			lineCounter++;
		}
        else
        {
            fprintf(stderr,"Unexpected type of kmers input file\n");
            return EXIT_FAILURE;
        }
	}
	
	manager->originalNumberOfKmers = kmersCounter;
	manager->maxNumberOfLeaves = manager->originalNumberOfKmers+1; //that's because we start counting from 1
    
	if(manager->includeReverseComplement)
		manager->maxNumberOfLeaves=2*manager->maxNumberOfLeaves;

    if(DEBUG_KWTREE) fprintf(stderr,"manager->maxNumberOfLeaves=%ld\n",(long)manager->maxNumberOfLeaves);

	return EXIT_SUCCESS;

    
}


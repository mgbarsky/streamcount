#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include "KeywordTree.h"

static const char *optString = "cp:m:k:d:f:h?";        

static const struct option longOpts[] = {
	{ "perform count", no_argument, NULL, 'c' },
	{ "pattern file name", required_argument, NULL, 'p' },
	{ "memory in MB", required_argument, NULL, 'm' },
	{ "input directory", required_argument, NULL, 'd' },
	{ "length of k-mer", required_argument, NULL, 'k' },
	{ "file of inputs", required_argument, NULL, 'f' },
	{ "help", no_argument, NULL, 'h' },
	{ NULL, no_argument, NULL, 0 }
};

void display_usage()
{
	printf("./streamcount [-c] -p pattern_file_name -k length of k-mers  [-m memory_in_MB] [-d input_directory] [-f file_with_input_file_names] [-h] list_of_input_file_names\n");
}

int main(int argc, char *argv[])
{
	
	int opt,i,c;
	char currentFileName[10000];
	char currentLine[10001];
	int longIndex=0;

	GlobalArgs globalArgs;
	/* Initialize globalArgs before we get to work. */
	globalArgs.patternFileName = NULL;  
    	globalArgs.memoryInMB = 1000;  //1 GB
	globalArgs.countOrNot =0; //default - not to count, but only to build the keyword tree from k-mers
	globalArgs.k =0; //has to specify
	globalArgs.isInputDirectory =0; //default - no directory name provided
	globalArgs.inputDirName =NULL; 
	globalArgs.isFileWithFileNames =0; //default - take file names from program arguments
	globalArgs.fileFileNames =NULL; 
	globalArgs.numInputFiles =0;
	globalArgs.inputFilesFromCmdLine =0;
	FILE *inputFP;

	opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
    	while( opt != -1 ) 
	{
        	switch( opt ) 
		{
            		case 'c':
                		globalArgs.countOrNot = 1; /* true */
                		break;
                
            		case 'p':
                		globalArgs.patternFileName = optarg;
                		break;
                
            		case 'm':
                		globalArgs.memoryInMB = atoi(optarg);
                		break;
                
            		case 'k':
                		globalArgs.k=atoi(optarg);
                		break;

			case 'd':
				globalArgs.isInputDirectory =1;  /* true */
                		globalArgs.inputDirName = optarg;
                		break;

			case 'f':
				globalArgs.isFileWithFileNames=1;  /* true */
                		globalArgs.fileFileNames = optarg;
                		break;
                
            		case 'h':   /* fall-through is intentional */
            		case '?':
                		display_usage();
                		break;  
                
            		default:
                		/* You won't actually get here. */
                		break;
        	}
        
        	opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
    	}

	if(globalArgs.countOrNot==1) //we are going to perform count - need input file names
 	{
		if(!globalArgs.isFileWithFileNames)  //if file with file names is not provided - file names are the rest of the command line arguments
		{
			globalArgs.inputFiles = argv + optind;
			globalArgs.numInputFiles = argc - optind;
			globalArgs.inputFilesFromCmdLine=1;
		}
		else  //we are going to read a file that contains input file names
		{
			
			sprintf(currentFileName,"%s", globalArgs.fileFileNames);
			
			if(!(inputFP= fopen ( currentFileName , "r" )))
			{
				printf("Could not open file \"%s\" for reading input file names \n", currentFileName);
				return 1;
			}
			
			
			while( fgets (currentLine, 10000, inputFP)!=NULL ) 
			{
				//avoid empty lines
				if(strlen(currentLine)>3)
					globalArgs.numInputFiles++;			
			}
			
			
			globalArgs.inputFiles = (char **)calloc(globalArgs.numInputFiles,sizeof (char *));
			rewind(	inputFP);
			i=0;
			while( fgets (currentLine, 10000, inputFP)!=NULL ) 
			{
				//avoid empty lines
				int lineLen = strlen(currentLine);
				if(lineLen>3)
				{
					globalArgs.inputFiles[i]=(char *)calloc(lineLen+1,sizeof (char ));
					//transfer file name from the line of input
					int endOfLine=0;
					for(c=0;c<lineLen && !endOfLine;c++)
					{
						if(currentLine[c]==10 ||currentLine[c]==32 ||currentLine[c]=='\0')
						{
							endOfLine=1;
							globalArgs.inputFiles[i][c]='\0';
						}
						else
						{
							globalArgs.inputFiles[i][c]=currentLine[c];
						}
					}
					if(!endOfLine)
						globalArgs.inputFiles[i][c]='\0';
					i++;	
				}		
			}
			fclose(inputFP);
		}
		if(globalArgs.numInputFiles==0)
		{
			printf("MANDATORY PARAMETER IS MISSING\n");
			printf("In order to count k-mers you need to specify a set of files where to count\n");
			display_usage();
			return 1;
		}		
	}

	if(globalArgs.patternFileName == NULL)
	{
		printf("MANDATORY PARAMETER IS MISSING\n");
		printf("You have to provide the name of a file with patterns to extract k-mers from it\n");
		display_usage();
		return 1;
	}

	if(!globalArgs.k)
	{
		printf("MANDATORY PARAMETER IS MISSING\n");
		printf("The value of k is not specified\n");
		display_usage();
		return 1;
	}
	
	
	if(performStreamCount(&globalArgs)!=0)
		return 1;

	//free memory
	if(!globalArgs.inputFilesFromCmdLine && globalArgs.inputFiles!=NULL)
	{
		for(i=0;i<globalArgs.numInputFiles;i++)
			free(globalArgs.inputFiles[i]); 

		free(globalArgs.inputFiles);
	}

	return 0;	
}

int performStreamCount(GlobalArgs* globalArgs)
{
	char currentFileName[MAX_PATH_LENGTH];
	int buildOrNot = 0;
	FILE *fp;
	
	//we are going to rebuild the keyword index if it does not exist
	sprintf(currentFileName, "%s_%d-mers_KWTREE_INFO", globalArgs->patternFileName, globalArgs->k);

	if(!(fp= fopen ( currentFileName , "rb" )))	
		buildOrNot=1;
	
	if(buildOrNot)
	{
		if(preprocessPatternsIntoKeywordTree(globalArgs)!=0)
		{
			printf("Error building keyword tree for the pattern set\n");
			return 1;
		}		
	}

	//if -c option is specified we are going to iterate over files in globalArgs.inputFiles and produce counts
	if(globalArgs->countOrNot==1)
	{
		if(countAll(globalArgs)!=0)
		{
			printf("Error counting\n");
			return 1;
		}
	}
	return 0;
}

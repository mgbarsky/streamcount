#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include "KeywordTree.h"

static const char *optString = "cp:m:k:d:f:o:h?ri:";        

static const struct option longOpts[] = {
	{ "perform count", no_argument, NULL, 'c' },
	{ "pattern file name", required_argument, NULL, 'p' },
	{ "memory in MB", required_argument, NULL, 'm' },
	{ "input directory", required_argument, NULL, 'd' },
	{ "output directory", required_argument, NULL, 'o' },
	{ "length of k-mer", required_argument, NULL, 'k' },
	{ "file of inputs", required_argument, NULL, 'f' },
	{ "help", no_argument, NULL, 'h' },
	{ "reverse complement", no_argument, NULL, 'r' },
	{ "input patterns type", required_argument, NULL, 'i' },
	{ NULL, no_argument, NULL, 0 }
};

void display_usage()
{
    static const char* USAGE_STRING =
	"./streamcount [options] -p pattern_file -k k_mer_length query_file1 query_file2 ..."
    "\n"
    "Options: \n"
    "        -h         Print this help and exit\n"
    "        -c         Count the number of occurrences of each pattern in the query files\n"
    "        -r         Include reverse-complement counts in the output\n"
    "        -f FILE    Read query file names from FILE\n"
    "        -i MODE    Patterns input MODE. MODE can be:\n"
    "                     0 - patterns are read from each line of the input\n"
    "                     1 - patterns are read from a FASTA file\n"
    "        -d DIR     Read input files from DIR\n"
    "        -o DIR     Write results to DIR\n"
    "        -m N       Use at most N megabytes for the patterns\n";
    printf("%s", USAGE_STRING);
}

int main(int argc, char *argv[])
{
	
	int opt,i,c;
	char currentFileName[10000];
	char currentLine[10001];
	int longIndex=0;

	GlobalArgs globalArgs;
	/* Initialize globalArgs with default values before we get to work. */
	globalArgs.patternFileName = NULL;  
    	globalArgs.memoryInMB = 1000;  //1 GB
	globalArgs.countOrNot =0; //default - not to count, but only to build the keyword tree from k-mers
	globalArgs.k =0; //has to specify
	globalArgs.isInputDirectory =0; //default - no directory name provided - use current directory
	globalArgs.inputDirName =NULL; 
	globalArgs.isFileWithFileNames =0; //default - take file names from program arguments (command line)
	globalArgs.fileFileNames =NULL; 
	globalArgs.numInputFiles =0;
	globalArgs.inputFilesFromCmdLine =0;
	globalArgs.isOutputDirectory =0; //no output directory provided - counts will be stored in a current directory of the input files
	globalArgs.includeReverseComplement=0; //need to specify -r if you want to include reverse complement
	globalArgs.inputType = INPUT_LINES;
	
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

			case 'o':
				globalArgs.isOutputDirectory =1;  /* true */
                		globalArgs.outputDirName = optarg;
                		break;

			case 'f':
				globalArgs.isFileWithFileNames=1;  /* true */
                		globalArgs.fileFileNames = optarg;
                		break;

                	case 'i':				
                		globalArgs.inputType = atoi(optarg);
                		break;

			case 'r':				
                		globalArgs.includeReverseComplement = 1;
                		break;

            		case 'h':   /* fall-through is intentional */
            		case '?':
                		display_usage();
                        exit(EXIT_SUCCESS);
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
			
			snprintf(currentFileName, MAX_PATH_LENGTH, "%s", globalArgs.fileFileNames);
			
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
			printf("In order to count k-mers you need to specify a set of files where to count\n\n");
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
	
	if(globalArgs.inputType<0 || globalArgs.inputType>2)
	{
		printf("INVALID TYPE OF PATTERNS INPUT\n");
		printf("The value should be one of: 0-lines, 1-file, 2-snip\n");
		display_usage();
		return 1;
	}
	
	if(process(&globalArgs)!=0)
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

int process(GlobalArgs* globalArgs)
{
	char currentFileName[MAX_PATH_LENGTH];
	int buildOrNot = 0;
	FILE *fp;
	
	//we are going to rebuild the keyword index if it does not exist
	snprintf(currentFileName, MAX_PATH_LENGTH, "%s_%d-mers_KWTREE_INFO", globalArgs->patternFileName, globalArgs->k);

	if(!(fp= fopen ( currentFileName , "rb" )))	
		buildOrNot=1;
	
	if(buildOrNot)
	{
		if(buildPatternIndex(globalArgs)!=0)
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

//Builds keyword tree with suffix links for each k-mer
//obtained from a set of patterns extracted from file globalArgs.patternFileName
//Converting each line of the input file (line length cannot exceed MAX_CHARS_PER_LINE) into k-mer substrings - depending on the type of the input
int 
buildPatternIndex(GlobalArgs *globalArgs)
{
	char patternFileName [MAX_PATH_LENGTH]; 
	KWTreeBuildingManager manager;
	int written;
	char currentFileName[MAX_PATH_LENGTH];
	FILE *outputFP;
	KWTreeInfo info[1];
	int64_t totalUniquePatterns=0;
	
	int memoryMB;	
	
	snprintf(patternFileName, MAX_PATH_LENGTH, "%s", globalArgs->patternFileName);	
	manager.k = globalArgs->k;
	manager.inputType = globalArgs->inputType;
	manager.includeReverseComplement = globalArgs->includeReverseComplement;	
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

	snprintf(currentFileName, MAX_PATH_LENGTH, "%s_%d-mers_KWTREE_INFO", patternFileName, manager.k);

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
	snprintf(currentFileName, MAX_PATH_LENGTH, "%s_%d-mers_KWTREE", patternFileName, manager.k);
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

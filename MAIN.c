#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include "General.h"
#include "StreamCount.h"
#include <unistd.h>
#include <getopt.h>

#define OPT_KMERS 1
#define OPT_KMERS_MULTILINE 2
#define OPT_INPUT_PLAIN_TEXT 3
#define OPT_NO_RC 4
#define OPT_PRINT_SEQ 5
#define OPT_REPEAT_MASK 6

static const char *shortOpts = "hk:i:m:e";        

static const struct option longOpts[] = {
	{ "kmers", required_argument, NULL, OPT_KMERS },
	{ "help", no_argument, NULL, 'h' },
	{ "k", required_argument, NULL, 'k' },
	{ "kmers-multiline", no_argument, NULL, OPT_KMERS_MULTILINE },
	{ "input", required_argument, NULL, 'i' },
	{ "input-plain-text", no_argument, NULL, OPT_INPUT_PLAIN_TEXT },
	{ "no-rc", no_argument, NULL, OPT_NO_RC },
	{ "mem", required_argument, NULL, 'm' },
	{ "printseq", no_argument, NULL, OPT_PRINT_SEQ },
	{ "repeat-mask-tofile", required_argument, NULL, OPT_REPEAT_MASK },
    { "e", no_argument, NULL, 'e' },
	{ NULL, no_argument, NULL, 0 }
};

void display_usage()
{
    static const char* USAGE_STRING =
    "Usage: ./streamcount [options] --kmers=<kmers_file> "
    "              --kmers=<kmers_file>                   extract k-mers from this file: required option\n"
    "Program for counting occurrences in one input file of all k-mers from a provided <kmers_file> - the input file can be provided as stdin \n"
    "\n"
    "Options: \n"
    "      -h,     --help                                      display this help and exit\n\n" 
    "  Input options: \n"  
    "      -k=<k>                                              length of each k-mer. \n"
    "                                                          Default: the length of the first valid line in <kmers_file> is used as <k>\n"
    "              --kmers-multiline                           extract k-mers from <kmers_file> treating it as one string.\n"
    "                                                          Default: k-mers are extracted from each separate line.\n"
    "      -i,     --input=<input_file>                        count k-mers in this file\n"
    "              --input-plain-text                          treat input as text lines. Default: FASTA-formatted input\n\n"
    "  Counting options: \n" 
    "              --no-rc                                     do not include count of reverse complement into final count of each k-mer. \n"
    "                                                          Default: include reverse complement.\n"
    "      -m,     --mem=<MEMORY_MB>                           amount of memory available for indexing k-mers. \n"
    "                                                          Used to estimate if can hold k-mers index before processing. Default: 4000MB \n\n"
    "  Output options: \n" 
    "              --printseq                                  print an original line of <kmers_file> before its count. \n"
    "                                                          Default: print counts only - each line of count(s) corresponds to the line in <kmers_file>\n"
    "              --repeat-mask-tofile=<repeat-mask-file>     for each k-mer prints to <repeat-mask-file> 0 or 1. 1 is printed if this k-mer is not unique (repeats) in the <kmers_file> \n"
    "      -e                                                  write a small file indicating the end of processing - to use with cluster. To use this option, input file name has to be specified. \n\n";
    
    printf("%s", USAGE_STRING);
}

//The program requires at least 2 input parameters
//1- name of file to extract k-mers from: by default each line is considered as a separate k-mer
//2 - name of file to count k-mers in: by default input type is FASTA, each sequence has its own header and is considered as a separate input string 
int main(int argc, char *argv[])
{
	char *kmersFileName=NULL;
    char *inputFileName=NULL;

    char *repeatsFileName=NULL;
   
    FILE *inputFP=NULL;
	FILE *kmersFP=NULL;
    FILE *repeatsFP=NULL;

    FILE *endFile=NULL;
    char endFileName[MAX_PATH_LENGTH];

	INT k=0; // default - the length of the first valid line
	int includeRC=1; //default - DNA with reverse complement
	int kmersInputType=KMERS_FROM_LINES; //default - extract k-mers from each line - no line crossing
    int fileInputType=INPUT_FASTA; //default - INPUT_FASTA, others: INPUT_LINES, INPUT_FILE
	int memoryMB=4000; //default
    int printKmersLine = 0;
    char currentLine[MAX_CHARS_PER_LINE];
    char validLine[MAX_CHARS_PER_LINE];
    KWTreeBuildingManager kwtreemanager; //bookkeper for indexing
    KWTCounterManager manager; //bookkeper for counting
   
    int indicateEnd = 0;

    INT *kmersCounts;
    int longIndex;
    long opt,i;

    opt = getopt_long( argc, argv, shortOpts, longOpts, &longIndex );
    while( opt != -1 ) 
	{
        switch( opt ) 
		{
            case OPT_KMERS:
                kmersFileName= optarg;
                break;                
           
            case 'k':
                k=atoi(optarg);
                break;

			case OPT_KMERS_MULTILINE:
				kmersInputType=KMERS_FROM_FILE;
                break;

			case 'i':
				inputFileName=optarg; 
                break;

			case OPT_INPUT_PLAIN_TEXT:
				fileInputType=INPUT_LINES;
                break;

            case OPT_NO_RC:				
                includeRC=0;
                break;

			case 'm':				
                memoryMB=atoi(optarg);
                break;
            
            case OPT_PRINT_SEQ:				
                printKmersLine = 1;
                break;

            case OPT_REPEAT_MASK:				
                repeatsFileName=optarg; 
                break;

            case 'e':
				if(inputFileName!=NULL)
                {
                    indicateEnd=1;
                    snprintf(endFileName, MAX_PATH_LENGTH, "%s_END", inputFileName);
                    if(!(endFile= fopen ( endFileName, "w" )))
                    {
                        fprintf(stderr,"Could not open file \"%s\" for writing end of program mark \n", endFileName);
                        return endProgram(EXIT_FAILURE, indicateEnd, endFile );
                    }
                } 
                break;

            case 'h':   /* fall-through is intentional */
            case '?':
                display_usage();
                return endProgram(EXIT_SUCCESS, 0, endFile );
                break;  
                
           default:
                /* You won't actually get here. */
                break;
        }
        
        opt = getopt_long( argc, argv, shortOpts, longOpts, &longIndex );
    }	
	
    //check that k-mers file is not NULL and try to open it
    if(kmersFileName == NULL)
    {
        fprintf(stderr,"Required parameter: --kmers name of file with k-mers \n");
        return endProgram(EXIT_FAILURE, indicateEnd, endFile );
    }

    if(!(kmersFP= fopen ( kmersFileName, "r" )))
    {
        fprintf(stderr,"Could not open file \"%s\" for extracting k-mers from it \n", kmersFileName);
        return endProgram(EXIT_FAILURE, indicateEnd, endFile );
    }

    //check that inputFileName is not NULL. Try to read from stdin
    if(inputFileName == NULL)
    {
        inputFP= stdin; 
    }
    else
    {
        if(!(inputFP= fopen ( inputFileName , "r" )))
        {
            fprintf(stderr,"Could not open file \"%s\" for counting k-mers in it \n", inputFileName);
            return endProgram(EXIT_FAILURE, indicateEnd, endFile );
        }
    }

    if(inputFP==NULL)
    {
        fprintf(stderr,"No input file provided: where am I supposed to count k-mers?\n");
        return endProgram(EXIT_FAILURE, indicateEnd, endFile );
    }

    //check if k==0 and if yes try to deduce from k-mers file
    if(k==0)
    {
        int OK=0;
        while( fgets (currentLine, MAX_CHARS_PER_LINE-10, kmersFP)!=NULL  && !OK) 
        {
            int lineLen = lenValidChars(currentLine,strlen(currentLine));
            if(lineLen >3)
            {
                OK=1;
                k=lineLen;
            }
        }

        if (k==0)
        {
            fprintf(stderr,"No valid lines in k-mers input file \"%s\" for extracting k-mers from it \n", kmersFileName);
            return endProgram(EXIT_FAILURE, indicateEnd, endFile );
        }
    }
    rewind(kmersFP);


    //now we have everything we need to perform counting
    //1. size of k
    //2. kmersFP pointing to the beginning of k-mers file
    //3. inputFP pointing to the input file
    //4. includeRC, kmersInputType, fileInputType, memoryMB
    //5. Output options: moreThanOneKmerInLine, printKmersLine, repeatsFileName

	//1. Setup fields of kwtreemanager - to pass reference to it rather than repeating arguments in each function
    kwtreemanager.k = k;
	kwtreemanager.inputType = kmersInputType;
	kwtreemanager.includeReverseComplement = includeRC;	
	
	kwtreemanager.maxSetSize = (memoryMB * 1000000)/sizeof(KWTNode); //how many KWT nodes can be held in this amount of memory
	
    //**************************************
    //TTT. first, extract k-mers and build a kw-tree - in RAM, without writing to disk
    //*************************************
    if(DEBUG_KMERS_EXTRACTION)         fprintf(stderr,"K-mers file is %s\n",kmersFileName);
    if(convertAllKmersIntoKWTreeReturnTree (kmersFP,  kmersInputType,  k,  includeRC,  memoryMB, &kwtreemanager)!=EXIT_SUCCESS)
		return endProgram(EXIT_FAILURE, indicateEnd, endFile );

    //free patterns array - dont need it anymore 
	for(i=0;i<kwtreemanager.estimatedNumberOfKmers;i++)
	{
		free(kwtreemanager.kmers[i]); 
	}
	free(kwtreemanager.kmers);	

    //**********************8
    //SSS. second, stream input file through existing KWTREE which is in kwtreemanager.KWtree
    //***************************
    //init fields of a new manager
    manager.numberOfKWTreeNodes = kwtreemanager.actualNumberOfKWTreeNodes;
	manager.k = k;
	manager.numberOfKWTreeLeaves =kwtreemanager.actualNumberOfLeaves;
    manager.KWTree = kwtreemanager.KWtree;
    manager.inputType = fileInputType;
    manager.inputFP = inputFP;

    //allocate memory for susbtring counters
	if(!(manager.substringCounts =(INT *)calloc(manager.numberOfKWTreeLeaves,sizeof (INT))))
	{
		fprintf(stderr,"Unable to allocate memory to hold %ld counters.\n",(long)manager.numberOfKWTreeLeaves );
		return endProgram(EXIT_FAILURE, indicateEnd, endFile );
	}
	
    fprintf(stderr,"________________________________________\n");
	fprintf(stderr,"Started counting k-mers in the input file  (may take some time)...\n");
    fprintf(stderr,"*********************************\n\n");
    if(streamAndCountOneFile(&manager)!=EXIT_SUCCESS)
			return endProgram(EXIT_FAILURE, indicateEnd, endFile );

if(PRINT_COUNTING)
{
    for(i=0; i< manager.numberOfKWTreeLeaves; i++)
    {   
        fprintf(stderr,"%ld\n",(long)manager.substringCounts[i]);
    }
}
    fprintf(stderr,"*************************************\n");
	fprintf(stderr,"Counting complete\n");
    fprintf(stderr,"________________________________________\n\n");

    //**************************************
    //GGG. convert substring counts into k-mer counts, using manager.substringCounts array and the mapping kmersInfo
    //*************************************	
    //allocate memory for k-mer counts
    if(DEBUG_COUNTING) fprintf(stderr,"Number of input k-mers =%ld\n",(long)kwtreemanager.originalNumberOfKmers);
    if(!(kmersCounts =(INT *)calloc(kwtreemanager.originalNumberOfKmers,sizeof (INT))))
	{
		fprintf(stderr,"Unable to allocate memory to hold %ld k-mer counters.\n",(long)kwtreemanager.originalNumberOfKmers );
		return endProgram(EXIT_FAILURE, indicateEnd, endFile );
	}      
        
    if(combineSubstringCountsIntoKmersCounts(manager.numberOfKWTreeLeaves,manager.substringCounts,
        kwtreemanager.originalNumberOfKmers, kwtreemanager.kmersInfo,   kmersCounts, includeRC )!=EXIT_SUCCESS)
    {
        fprintf(stderr,"Unexpected error producing final k-mer counts.\n" );
		return endProgram(EXIT_FAILURE, indicateEnd, endFile );
    }

    //**************************************
    //OOO. Output
    //*************************************	
       
    if(printKmersLine)
    {        
        rewind(kmersFP);
        //read first valid line
        if(nextValidLineTextFile (kmersFP, k, validLine)!=EXIT_SUCCESS)        
        {
            fprintf(stderr,"Unexpected error: no valid line found for the 0 line of k-mers.\n" );
		    return endProgram(EXIT_FAILURE, indicateEnd, endFile );
        }
    }

    //now print counts to stdout -
    int currentLineNumber=0;
    if(DEBUG_COUNTING) fprintf(stderr,"N of k-mers = %ld\n",(long)kwtreemanager.originalNumberOfKmers);
    for(i=0; i< kwtreemanager.originalNumberOfKmers; i++)
    {
        if(currentLineNumber != kwtreemanager.kmersInfo[i].lineNumber)  //next line
        {
            fprintf(stdout,"\n");
            currentLineNumber=kwtreemanager.kmersInfo[i].lineNumber;

            if(printKmersLine)
            {
                //read the next line - if there is the next line
                if(i< kwtreemanager.originalNumberOfKmers-1)
                {
                    if(nextValidLineTextFile (kmersFP, k, validLine)!=EXIT_SUCCESS)        
                    {
                        fprintf(stderr,"Unexpected error: no valid line found for the %d-th line of k-mers.\n",currentLineNumber );
		                return endProgram(EXIT_FAILURE, indicateEnd, endFile );
                    }
                }
            }
        }


        INT absCount = (kmersCounts[i]<0)?- kmersCounts[i]: kmersCounts[i];
        
        if(kwtreemanager.kmersInfo[i].startPosInLine == 0)  //print first in line
        {
            if(printKmersLine)
            {
                fprintf(stdout,"%s,%ld",validLine,(long)absCount);                
            }
            else
            {
                fprintf(stdout,"%ld",(long)absCount);
            }  
        }
        else
        {            
             fprintf(stdout,",%ld",(long)absCount);            
        }       
    }
    
    fclose(kmersFP);
    //now print repeats masks if required in the same format as counts
    //open repeats file for output
    if(repeatsFileName!=NULL) //output repeat mask
    {
        if(!(repeatsFP= fopen ( repeatsFileName, "w" )))
        {
            fprintf(stderr,"Could not open file \"%s\" for writing repeats masks for input k-mers \n", repeatsFileName);
            return endProgram(EXIT_FAILURE, indicateEnd, endFile );
        }

        int currentLineNumber=0;
        for(i=0; i< kwtreemanager.originalNumberOfKmers; i++)
        {
            if(currentLineNumber != kwtreemanager.kmersInfo[i].lineNumber)  //next line
            {
                fprintf(repeatsFP,"\n");
                currentLineNumber=kwtreemanager.kmersInfo[i].lineNumber;
            }
            int repeating=(kmersCounts[i]<0)?1:0;            
            
            if(kwtreemanager.kmersInfo[i].startPosInLine == 0)  //print first in line
            {
                fprintf(repeatsFP,"%d",repeating);
            }
            else
            {
                fprintf(repeatsFP,",%d",repeating);
            }           
        }
    }

    //free all arrays in memory	
	free(kwtreemanager.kmersInfo);

	//free KWtree
	if(kwtreemanager.KWtree)
		free(kwtreemanager.KWtree);
    //free counts memory
	free(manager.substringCounts);	
    free(kmersCounts);

    
	return endProgram(EXIT_SUCCESS, indicateEnd, endFile );
}

int endProgram(int exitStatus, int indicateEnd, FILE *endFile )
{
    if(indicateEnd)
    {
        fprintf(endFile,"%d",exitStatus);
        fclose(endFile);
    }
    return exitStatus;
}

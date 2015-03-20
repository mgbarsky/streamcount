/*
This will extract a predefined number (about) of random valid DNA k-mers from a fasta file or from a file containing lines of text
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include "streamcount.h"

#include <zlib.h>
#include "kseq.h"

#define DEBUG 1

// Initialize KSEQ reader
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	SC_INT k;
	SC_INT howMany;
	char *inputFileName=NULL;
	
	FILE * inputFP;
	char currentLine[MAX_CHARS_PER_LINE];
    
    int inputType=INPUT_FASTA; //default
	
	long totalLines;
	//long kmersInOneLine;
	long totalKmersInFile;

	//double frequency;
   // double kmersPerLine;
	char ** kmers;

	SC_INT totalExtracted;
  
    gzFile gzFP; 
    kseq_t* seq;

    SC_INT i;
    
	if(argc<5)
	{
		fprintf(stderr,"./sc_extractrandomkmers <size of k-mers> <howMany> <input-plain-text:0-no, 1-yes> <input file name>  \n");
		return 1;
	}

	k=atoi(argv[1]);
	howMany = atoi(argv[2]);
    inputType=atoi(argv[3]);    
    inputFileName=argv[4];
    
    if(DEBUG)
    {
        fprintf(stderr, "Extracting %ld k-mers of length %ld from file '%s' \n",(long)howMany, (long)k,inputFileName);
    }
        
    if(!(inputFP= fopen ( inputFileName , "r" )))
    {
        fprintf(stderr,"Could not open file \"%s\" for estimating number of k-mers in file \n", inputFileName);
        return EXIT_FAILURE;
    }

    //allocate memory to hold a set of extracted k-mers
	//they should fit in RAM - subsequent processing relies on that
	if(!(kmers =(char **)calloc(howMany,sizeof (char *)))) 
	{
		fprintf(stderr,"Unable to allocate memory for kmers array to be extracted from a file.\n");
		return EXIT_FAILURE;
	}

	for( i=0;i<howMany;i++)
	{
		if(!(kmers[i] =(char *)calloc(k+1,sizeof (char *))))
		{
			fprintf(stderr,"Unable to allocate memory for kmer %ld.\n",(long)i);
			return EXIT_FAILURE;
		}
	}

	//1. estimate how many total k-mers can be extracted from the input file
    //read file line-by-line and add up total potential k-mers, as well as total number of lines
    totalKmersInFile=0;
    totalLines =0;

    if(inputType == INPUT_FASTA)
    { 
        if(!( gzFP = gzdopen ( fileno(inputFP) , "r" )))
	    {
		    fprintf(stderr,"Could not open input file as gz for reading\n");
		    return EXIT_FAILURE;
	    }

        seq = kseq_init(gzFP);

        // read sequences
        while(kseq_read(seq) >= 0 )
        {            
            int lineLen = lenValidChars(seq->seq.s, seq->seq.l);
            int kmersInLine = MAX(0, lineLen - k +1);	

            totalKmersInFile+= kmersInLine;
            totalLines++;
        }

        kseq_destroy(seq);
	    gzclose(gzFP);
    }
    else
    {
        while( fgets (currentLine, MAX_CHARS_PER_LINE-10, inputFP)!=NULL) 
	    {
		    int lineLen = lenValidChars(currentLine, strlen(currentLine));		
            
            int kmersInLine = MAX(0, lineLen - k +1);	

            totalKmersInFile+= kmersInLine;
            totalLines++;
	    }
    }

	if(DEBUG) fprintf(stderr,"Total potential k-mers = %ld, total lines=%ld\n", (long)totalKmersInFile, (long)totalLines);
    if(totalKmersInFile<=0)
    {
        fprintf(stderr,"No valid DNA k-mers found in the input file. Check your input file (and input type: 0 -fasta, 1 -text) which should contain lines of valid DNA characters\n");
        return EXIT_FAILURE;
    }

	
    fclose(inputFP);
    if(!(inputFP= fopen ( inputFileName , "r" )))
    {
        fprintf(stderr,"Could not open file \"%s\" for extracting rand kmers \n", inputFileName);
        return EXIT_FAILURE;
    }
    
	
	totalExtracted = howMany;

//EXTRACT k-mers
	if(extractKmers(inputFP, kmers,  k, &totalExtracted, totalLines, inputType)!=EXIT_SUCCESS)
	{
		fprintf(stderr,"Error occurred in 'fillKMersArray' \n");
		return EXIT_FAILURE;
	}
	fclose(inputFP);

	fprintf(stderr,"Extracted a set of %ld random k-mers from file %s\n",(long)totalExtracted,inputFileName);

	//write it to a stdout - check that all characters are DNA characters    
    SC_INT totalValid=0;
	for(i=0; i<totalExtracted; i++)
	{
        //check that all chars are valid DNA characters
        int lineLen = lenValidChars(kmers[i],strlen(kmers[i]));
        if(lineLen == k)
        {
		    fprintf(stdout,"%s\n", kmers[i]);
            totalValid++;
        }		
	}

    fprintf(stderr,"Written to stdout a set of %ld random k-mers with valid DNA alphabet\n",(long)totalValid);

	//free kmers array 
	for(i=0;i<howMany;i++)
	{
		free(kmers[i]); 
	}
	free(kmers);	
	return 0;
}

int extractKmers(FILE *inputFP, char **kmers, SC_INT k, SC_INT *totalKMers, SC_INT totalLines, int inputType)
{
	char currentLine[MAX_CHARS_PER_LINE];
	int lineLen;
	
	int counterTotalExtracted=0;
    SC_INT extractMaxNumber = *totalKMers;

    //now we want to know how many to extract from each line
    double howManyFromEachLine = (double)extractMaxNumber/totalLines;

    SC_INT skipLines = 0;
    SC_INT randFromLine = 1;

    if(howManyFromEachLine < 1)  //evenly distribute skipping of lines - more lines that the number of k-mers
    {
        skipLines = (SC_INT) (totalLines/extractMaxNumber +0.5);
    }
    else if(howManyFromEachLine > 1) //more than 1 k-mer from each line
    {
        randFromLine = (SC_INT) (howManyFromEachLine+0.5);
    }

    long skippedLinesCounter = 0;
    long extractedFromLineKmersCounter = 0;
	
    if(DEBUG) fprintf(stderr, "Want %ld k-mers from %ld lines, skip each %ld lines, extract %ld k-mer(s) from each line\n",
                            (long)extractMaxNumber, (long)totalLines, (long)skipLines, (long)randFromLine); 
    gzFile gzFP; 
    kseq_t* seq;

		
    //now it depends on the type of the input file
    if(inputType == INPUT_FASTA)
    { 
        if(!( gzFP = gzdopen ( fileno(inputFP) , "r" )))
	    {
		    fprintf(stderr,"Could not open input file as gz for the second time for reading\n");
		    return EXIT_FAILURE;
	    }

        seq = kseq_init(gzFP);

        // read sequences
        while(kseq_read(seq) >= 0 && counterTotalExtracted <extractMaxNumber)
	    {
		    if(seq->seq.l >=k && (!skipLines || skippedLinesCounter>=skipLines))  //use this line
            {
                int maxStartPos = seq->seq.l-k+1;
                
                if(randFromLine == 1)  //only 1 random k-mer is required in each line
                {
                    int randStart =  rand() % (maxStartPos+1);         // randStart in the range 0 to maxStartPos

                    memcpy(kmers[counterTotalExtracted],&seq->seq.s[randStart],k);
					kmers[counterTotalExtracted][k]='\0';				
					counterTotalExtracted++;				
                }
                else  //continue until extracted randFromLine k-mers from each line
                {
                    while (extractedFromLineKmersCounter < randFromLine && counterTotalExtracted <extractMaxNumber)
                    {
                        int randStart =  rand() % maxStartPos;         // randStart in the range 0 to maxStartPos

                        memcpy(kmers[counterTotalExtracted],&seq->seq.s[randStart],k);
					    kmers[counterTotalExtracted][k]='\0';				
					    counterTotalExtracted++;

                        extractedFromLineKmersCounter++;	
                    }

                    extractedFromLineKmersCounter=0;
                }   

                skippedLinesCounter=0; //reset                                           
            }  
            else //skip this line 
            {
                skippedLinesCounter++;
            }        	
        }

        kseq_destroy(seq);
	    gzclose(gzFP);       
    }
    else if (inputType == INPUT_LINES)
    {
	    while( fgets (currentLine, MAX_CHARS_PER_LINE-10, inputFP)!=NULL && counterTotalExtracted <extractMaxNumber) 
	    {
		    lineLen = lenValidChars(currentLine,strlen(currentLine));

            if(lineLen >=k && (!skipLines || skippedLinesCounter>=skipLines))  //use this line
            {
                int maxStartPos =lineLen-k+1;
                
                if(randFromLine == 1)  //only 1 random k-mer is required in each line
                {
                    int randStart =  rand() % (maxStartPos+1);         // randStart in the range 0 to maxStartPos

                    memcpy(kmers[counterTotalExtracted],&currentLine[randStart],k);
					kmers[counterTotalExtracted][k]='\0';				
					counterTotalExtracted++;				
                }
                else  //continue until extracted randFromLine k-mers from each line
                {
                    while (extractedFromLineKmersCounter < randFromLine && counterTotalExtracted <extractMaxNumber)
                    {
                        int randStart =  rand() % maxStartPos;         // randStart in the range 0 to maxStartPos

                        memcpy(kmers[counterTotalExtracted],&currentLine[randStart],k);
					    kmers[counterTotalExtracted][k]='\0';				
					    counterTotalExtracted++;

                        extractedFromLineKmersCounter++;	
                    }

                    extractedFromLineKmersCounter=0;
                }   

                skippedLinesCounter=0; //reset                                           
            }  
            else //skip this line 
            {
                skippedLinesCounter++;
            }  
        }
	}    
    else if(inputType == INPUT_FILE)
    {
        fprintf(stderr,"To be done later\n"); //TBD
        return EXIT_FAILURE;
    }	
    else
    {
        fprintf(stderr,"Invalid input type: %d\n",inputType);
        return EXIT_FAILURE;
    }   	
	
	*totalKMers = counterTotalExtracted;
	return EXIT_SUCCESS;
}


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <zlib.h>
#include "kseq.h"
#include "KeywordTree.h"

// Initialize KSEQ reader
KSEQ_INIT(gzFile, gzread)

//second part of stream and count program
//assumes that the input pattern set has been preprocessed into a keyword tree and serialized to a file <%s_%d-mers_KWTREE_INFO>
//parameters -  names of input files stored in globalArgs.inputFiles, total number of input files stored in globalArgs.numInputFiles
int countAll(GlobalArgs *globalArgs)
{
	KWTCounterManager manager; //bookkeper
	char currentFileName[MAX_PATH_LENGTH];
	KWTreeInfo newinfo[1];
	FILE *kwInputFP;
    FILE *outputFP;
	int read, written;
	int i;
	int f;	

	//try to read KWtree info from file
	snprintf(currentFileName, MAX_PATH_LENGTH, "%s_%d-mers_KWTREE_INFO", globalArgs->patternFileName, globalArgs->k);
	if(!(kwInputFP = fopen ( currentFileName , "rb" )))
	{
		printf("Could not open file \"%s\" for reading saved KWtree info \n", currentFileName);
		return 1;
	}		
	
	read = fread (newinfo,sizeof(KWTreeInfo),1,kwInputFP);
	if(read!=1)
	{
		printf("Error reading KWTinfo from file\n");
		fclose(kwInputFP);
		return 1;
	}
	fclose(kwInputFP);

	//print new tree parameters read from file
	printf("Read from file %s the following: totalNodes = %d, k=%d, totalPatterns=%d\n",
		currentFileName,newinfo[0].treeSlotsNum,newinfo[0].k,newinfo[0].totalLeaves);

	//init general fields of the manager
	manager.treeSlotsNum = newinfo[0].treeSlotsNum;
	manager.k = newinfo[0].k;
	manager.totalPatterns = newinfo[0].totalLeaves;
	
	//allocate memory for KWT tree
	if(!(manager.KWTree = (KWTNode *)calloc(manager.treeSlotsNum,sizeof (KWTNode))))
	{
		printf("Unable to allocate memory to load keyword tree.\n");
		return 1;
	}
	
	//open file for reading a tree into RAM
	snprintf(currentFileName, MAX_PATH_LENGTH, "%s_%d-mers_KWTREE", globalArgs->patternFileName, manager.k);
	if(!(kwInputFP = fopen ( currentFileName , "rb" )))
	{
		printf("Could not open file \"%s\" for reading saved KWtree \n", currentFileName);
		return 1;
	}

	//read tree from file
	read= fread (&(manager.KWTree[0]), sizeof(KWTNode), manager.treeSlotsNum, kwInputFP);
	if(read!=manager.treeSlotsNum)
	{
		printf("Error reading KWT tree from file: wanted to read %d nodes but read %d\n", manager.treeSlotsNum, read);
		return 1;
	}
	fclose(kwInputFP);

	//allocate memory for counters
	if(!(manager.patternCounts =(UINT *)calloc(manager.totalPatterns,sizeof (UINT))))
	{
		printf("Unable to allocate memory to hold %d counters.\n",manager.totalPatterns );
		return 1;
	}
	
	for(f=0;f<globalArgs->numInputFiles;f++)
	{
		//counting in a specific input file 
		if(globalArgs->isInputDirectory)
			snprintf(manager.inputFileName, MAX_PATH_LENGTH, "%s//%s",globalArgs->inputDirName,globalArgs->inputFiles[f]);
		else
			snprintf(manager.inputFileName, MAX_PATH_LENGTH, "%s",globalArgs->inputFiles[f]);

		//reset counters from the previous file
		if(f>0)
		{
			for(i=0;i<manager.totalPatterns;i++)
				manager.patternCounts[i]=0;
		}

		//now we are ready - KWtree is in memory
		//stream each line of the input file through kwtree and update counts
		if(streamAndCountOneFile(&manager)!=0)
			return 1;	

		//output counters to a file		
		if(globalArgs->isOutputDirectory)
			snprintf(currentFileName, MAX_PATH_LENGTH, "%s//%s_%d-mers_COUNTS", globalArgs->outputDirName, globalArgs->inputFiles[f],manager.k);
		else
			snprintf(currentFileName, MAX_PATH_LENGTH, "%s_%d-mers_COUNTS", manager.inputFileName, manager.k);
		if(!(outputFP= fopen ( currentFileName , "wb" )))
		{
			printf("Could not open file \"%s\" for writing counts \n", currentFileName);
			return 1;
		}

		written = fwrite (&(manager.patternCounts[0]) , sizeof(UINT), manager.totalPatterns, outputFP);

		if(written!=manager.totalPatterns)
		{
			printf("Not all counters have been written: wanted to write %d, and wrote %d \n",manager.totalPatterns,written);
			return 1;
		}
		fclose(outputFP);
		
	}

	//free allocated memory
	free(manager.patternCounts);
	free(manager.KWTree);
	printf("All files have been processed\n");
	return 0;
}

int streamAndCountOneFile(KWTCounterManager *manager)
{
	gzFile inputFP;	
    kseq_t* seq;

	//open file to read lines
	if(!( inputFP = gzopen ( manager->inputFileName , "r" )))
	{
		printf("Could not open input file \"%s\" for reading\n", manager->inputFileName);
		return 1;
	}

    // initialize reader
    seq = kseq_init(inputFP);

    // read sequences
    while(kseq_read(seq) >= 0)
	{
		if(streamOneStringUnchanged(manager, seq->seq.s, seq->seq.l) != 0)
		{
			gzclose(inputFP);
			return 1;
		}	
    }

    kseq_destroy(seq);
	gzclose(inputFP);
	return 0;	
}

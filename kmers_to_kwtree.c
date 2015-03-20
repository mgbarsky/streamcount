/*
This will extract all possible k-mers from the input file
and convert them into an in-mem keyword tree
saving information about k-mers - their mapping to input lines, and to the corresponding reverse complement
*/
#include "kmers_to_kwtree.h"

int convertAllKmersIntoKWTreeReturnTree (FILE *kmersFP, int inputType, 
                                        SC_INT k, int includeRC, SC_INT memoryMB, 
                                        KWTreeBuildingManager *manager) {
	SC_INT estimatedNumberOfKmers=0;
	SC_INT estimatedTotalNumberOfKmers; //twice estimatedNumberOfKmers if include rc is true	
	SC_INT i;
	
   	//this will read input file and estimate how many k-mers it is possible to extract - return into estimatedNumberOfKmers
	if(collectKmerInputStats(kmersFP,manager->k,manager->inputType, &estimatedNumberOfKmers) != EXIT_SUCCESS)	{
		fclose(kmersFP);
		return EXIT_FAILURE;
	}

	if(DEBUG_KMERS_EXTRACTION) fprintf(stderr,"Estimated total %ld non-rc k-mers of size %ld each\n", (long)estimatedNumberOfKmers,(long)(manager->k)); 
	
	if (estimatedNumberOfKmers == 0)	{
		fprintf(stderr,"No valid characters found in the pattern input file, or the k-mer size exceeds the line length\n");
		return EXIT_FAILURE;
	}	
	
	//check if we have sufficient memory to hold KWtree
	estimatedTotalNumberOfKmers = (manager->includeReverseComplement) ? 2*estimatedNumberOfKmers : estimatedNumberOfKmers;
	manager->estimatedNumberOfLeaves = estimatedTotalNumberOfKmers+1;

	if(estimatedTotalNumberOfKmers*(manager->k+1) < manager->maxSetSize && estimatedTotalNumberOfKmers*(manager->k+1) < MAX_SC_INT)
		manager->estimatedNumberOfKWTreeNodes = estimatedTotalNumberOfKmers*(manager->k+1);
	else	{
		fprintf(stderr,"Pattern set can not be preprocessed in the amount of the available memory:\n");
		fprintf(stderr,"We can hold maximum %ld keyword tree nodes, but the tree may require %ld nodes\n",
				(long)manager->maxSetSize, (long)estimatedTotalNumberOfKmers*(manager->k+1));
		return EXIT_FAILURE;
	}
	
	//allocate an array to hold all kmers and their info
	if ( !(manager->kmers = (char **)calloc(estimatedNumberOfKmers,sizeof (char *))) )	{
		fprintf(stderr,"Unable to allocate memory for kmers array.\n");
		return EXIT_FAILURE;
	}

	for (i=0; i<estimatedNumberOfKmers; i++)	{
		if( !(manager->kmers[i] = (char *)calloc(manager->k+1,sizeof (char *))) )  {
			fprintf(stderr,"Unable to allocate memory for kmer %ld.\n",(long)i);
			return EXIT_FAILURE;
		}
	}
	
	if( !(manager->kmersInfo = (KmerInfo *)calloc(estimatedNumberOfKmers,sizeof (KmerInfo))) ) 	{
		fprintf(stderr,"Unable to allocate memory for kmers info array.\n");
		return EXIT_FAILURE;
	}

    //remember estimatedNumberOfKmers to free these arrays later
    manager->estimatedNumberOfKmers = estimatedNumberOfKmers;

	rewind(kmersFP); 

	//after this, the actual number of valid k-mers to process is in manager->originalNumberOfKmers
	if(fillKmersArrayAndInfo(manager, kmersFP, manager->kmers, manager->kmersInfo, estimatedNumberOfKmers)!=EXIT_SUCCESS){
		fprintf(stderr,"Error generating kmers array\n");
		fclose(kmersFP);
		return EXIT_FAILURE;
	}

	if (DEBUG_KMERS_EXTRACTION) fprintf(stderr,"Preprocessing of the input file into an array of kmers complete\n");

	//allocate memory for KWtree
	SC_INT maxNumberOfNodes = (manager->maxNumberOfLeaves) * (manager->k) + 1; //each of k chars can be potentially a node plus a root node
	if(!(manager->KWtree = (KWTNode *) calloc (maxNumberOfNodes,sizeof (KWTNode))) ){
		fprintf(stderr,"Unable to allocate memory for keyword tree.\n");
		return EXIT_FAILURE;
	}
    
	if (DEBUG_KWTREE)    {
        unsigned long allMemInBytes=maxNumberOfNodes*sizeof (KWTNode);
        fprintf(stderr,"Allocated memory %lu to hold %ld KW tree nodes with max %ld number of leaves \n", 
                allMemInBytes,(long) maxNumberOfNodes,(long)manager->maxNumberOfLeaves);
    }

	//pre-process patterns into a keyword tree - there duplicate k-mers get removed by design
	if (buildKeywordTree (manager, manager->kmers, manager->kmersInfo) != EXIT_SUCCESS)
		return EXIT_FAILURE;
    
    if (manager->saveTree)
    {
        saveTreeAndMapping (manager->kwtreeFileName, manager->KWtree, manager->actualNumberOfKWTreeNodes, manager->kmersInfo, manager->originalNumberOfKmers);
    }
	if(DEBUG_KWTREE) 
        fprintf(stderr,"BuildKeywordTree complete. totalNodes = %ld, k=%ld, totalLeaves=%ld\n", 
            (long)manager->actualNumberOfKWTreeNodes,(long)manager->k,(long)manager->actualNumberOfLeaves);	
	
    fprintf(stderr,"Preprocessed %ld-mers into index of %ld unique strings.\n",  (long)manager->k, (long)manager->actualNumberOfLeaves );   
	return EXIT_SUCCESS;
}

/* This serves to count total number of possible k-mers in the input
In order to allocate constant memory. Result is returned in *estimatedNumberOfKmers
*/
int collectKmerInputStats(FILE *inputFP, SC_INT k, int inputType, SC_INT *estimatedNumberOfKmers){
	char currentLine[MAX_CHARS_PER_LINE];
	SC_INT kmersCount=0;
	
	while( fgets (currentLine, MAX_CHARS_PER_LINE-10, inputFP)!=NULL ) 	{
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

/*this reads the input file and fills-in a pre-allocated k-mer array with sequences
as well as stores information about each k-mer such as 
    line number, 
    start position in the corresponding line 
*/
int fillKmersArrayAndInfo(KWTreeBuildingManager* manager, FILE *inputFP, 
	char **kmers, KmerInfo *kmersInfo, SC_INT maxPossibleNumberOfKmers)
{
	char currentLine[MAX_CHARS_PER_LINE];
	
	int i;
	int kmersCounter=0;
	int lineCounter = 0;
	
    if(manager->inputType == KMERS_FROM_LINES) {  //k-mers are extracted from each separate line of the input file
        while( fgets (currentLine, MAX_CHARS_PER_LINE-10, inputFP)!=NULL ) 	{
		
			int lineLen = lenValidChars(currentLine,strlen(currentLine));
			if(lineLen >= manager->k)	{  //only works with input file where each line contains at least k consecutive valid chars
				for(i=0;i<lineLen-manager->k+1;i++)	{
					memcpy(kmers[kmersCounter],&currentLine[i],manager->k);
					kmers[kmersCounter][manager->k]='\0';
					kmersInfo[kmersCounter].lineNumber = lineCounter;
					kmersInfo[kmersCounter].startPosInLine = i;
					kmersInfo[kmersCounter].repeated = 0; //initialize repetitions flag
					kmersCounter++;

					if(kmersCounter > maxPossibleNumberOfKmers) { //to avoid memory overflow
						fprintf(stderr,"Unexpected error while extracting k-mers: too many k-mers found !!\n");
						return EXIT_FAILURE;
					}
				}				
			}
            else  {
                fprintf(stderr,"Line %d contains non-DNA character or is less than k=%ld: %s\n",lineCounter,(long)manager->k,currentLine);
                return EXIT_FAILURE;
            }
			lineCounter++;			
		}
    }
	else if (manager->inputType == KMERS_FROM_FILE)	{  //treats input as one big string
    //in this case we are going to read input into a concatenated buffer, and record how many characters in each line into another buffer of line lengths
        fseek(inputFP, 0, SEEK_END);
        SC_INT totalChars = 0;
        SC_INT maxChars = (SC_INT)ftell(inputFP);
        rewind (inputFP);
    
        char *concatenatedInput = NULL;
        if(!(concatenatedInput = (char *)calloc(maxChars,sizeof (char))))	{
		    fprintf(stderr,"Unable to allocate memory to buffer %ld chars of the input for k-mer extraction.\n",(long)maxChars );
		    return EXIT_FAILURE;
	    }  

        SC_INT *inputLineLengths = NULL;
        if(!(inputLineLengths = (SC_INT *)calloc(maxChars,sizeof (SC_INT))))	{
		    fprintf(stderr,"Unable to allocate memory to buffer %ld line lengths in k-mer extraction.\n",(long)maxChars );
		    return EXIT_FAILURE;
	    }     
    
        //read the entire input into the allocated buffer, keep track to line numbers
        while( fgets (currentLine, MAX_CHARS_PER_LINE-10, inputFP)!=NULL ) 	{
		
			int lineLen = lenValidChars (currentLine, strlen(currentLine));
            inputLineLengths [ lineCounter++] = lineLen;
            
            //append current line to the concatenated string in buffer
            memcpy (&concatenatedInput[totalChars], currentLine, lineLen);
            totalChars += lineLen;
        }

        //start extracting k-mers from concatenated input
        //they can start anywhere except the last k-1 positions
        lineCounter = 0;
        SC_INT startPosInLine = 0;
        SC_INT currentLineLength = inputLineLengths [lineCounter];
        
        for (i=0; i < totalChars - (manager->k -1); i++) {
            memcpy (kmers[kmersCounter], &concatenatedInput[i], manager->k);
		    kmers[kmersCounter][manager->k] = '\0'; //terminate string

            kmersInfo[kmersCounter].lineNumber = lineCounter;  //starts in the previous line
			kmersInfo[kmersCounter].startPosInLine = startPosInLine++;

			kmersInfo[kmersCounter].repeated = 0; //initialize repetitions flag
            
            kmersCounter++;

            if (startPosInLine == currentLineLength ) {
                startPosInLine = 0;
                lineCounter ++;
                currentLineLength = inputLineLengths [lineCounter];
            }
        } 

        if (concatenatedInput) free (concatenatedInput);
        if (inputLineLengths ) free (inputLineLengths);
	}
	else  {
        fprintf(stderr,"Unexpected type of kmers input file\n");
        return EXIT_FAILURE;
    }

	manager->originalNumberOfKmers = kmersCounter;
    
    fprintf(stderr,"Extracted total %ld %ld-mers from k-mers file\n", 
            (long)kmersCounter,(long)manager->k);
    
	manager->maxNumberOfLeaves = manager->originalNumberOfKmers+1; //that's because we start counting from 1
    
	if(manager->includeReverseComplement){
		manager->maxNumberOfLeaves = 2 * manager->maxNumberOfLeaves;
        fprintf(stderr,"Count for each k-mer will contain count of its reverse complement.\n"); 
    }
    if (DEBUG_KMERS_EXTRACTION) {
        for (int i=0; i< 40; i++)
            fprintf(stderr,"%s\t%ld\t%ld\n", kmers[i], (long)kmersInfo[i].lineNumber, (long)kmersInfo[i].startPosInLine);
    }
    if(DEBUG_KWTREE) fprintf(stderr,"manager->maxNumberOfLeaves=%ld\n", (long)manager->maxNumberOfLeaves);

	return EXIT_SUCCESS;
}

int  saveTreeAndMapping (char *kwtreeFileName, KWTNode *KWtree, SC_INT totalNodes, KmerInfo *kmersInfo, SC_INT totalKmers)
{
    int i=0;
    
    char infoFileName [MAX_PATH_LENGTH];
    snprintf ( infoFileName, MAX_PATH_LENGTH, "%s_KWTREE_INFO.csv", kwtreeFileName);
    
    FILE *fpinfo;
    fpinfo=fopen(infoFileName, "w");
    if(!(fpinfo= fopen ( infoFileName, "w" ))) {
        fprintf(stderr,"Could not open file \"%s\" for writing kwtree info \n", infoFileName);
        return EXIT_FAILURE;
    }
    
    for( i = 0; i < totalKmers; i++) {
        fprintf(fpinfo, "%d,%d,%d,%d\n",(int)(kmersInfo[i].counterID),(int)(kmersInfo[i].rcCounterID),(int)(kmersInfo[i].startPosInLine),(int)(kmersInfo[i].lineNumber));
    }
    fclose (fpinfo);


    FILE *fp;
    if(!(fp= fopen ( kwtreeFileName, "wb" ))) {
        fprintf(stderr,"Could not open file \"%s\" for writing kwtree \n", kwtreeFileName);
        return EXIT_FAILURE;
    }

    size_t written = fwrite(KWtree, sizeof(KWTNode), totalNodes, fp);
    if (written != (size_t)totalNodes) {
        fprintf(stderr,"Not all k-mers have been written \n");
        return EXIT_FAILURE;
    }
    fclose (fp);

    return EXIT_SUCCESS;
}





















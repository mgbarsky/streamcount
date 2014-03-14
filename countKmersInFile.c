#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

#include "General.h"
#include "StreamCount.h"

#include <zlib.h>
#include "kseq.h"

// Initialize KSEQ reader
KSEQ_INIT(gzFile, gzread)

int streamAndCountOneFile(KWTCounterManager *manager)
{	
	
	char currentLine[MAX_CHARS_PER_LINE];
	gzFile gzFP; 
    kseq_t* seq;
    INT validLines=0;

	//now it depends on the type of the input file
    if(manager->inputType == INPUT_FASTA)
    { 
        if(!( gzFP = gzdopen ( fileno(manager->inputFP) , "r" )))
	    {
		    fprintf(stderr,"Could not open input file as gz for reading\n");
		    return EXIT_FAILURE;
	    }

        seq = kseq_init(gzFP);

        // read sequences
        while(kseq_read(seq) >= 0)
	    {
		    if(streamOneString(manager->KWTree, seq->seq.s, seq->seq.l,&manager->substringCounts[0] )!=EXIT_SUCCESS )
		    {
			    gzclose(gzFP);
			    return EXIT_FAILURE;
		    }
            validLines++;	
        }

        kseq_destroy(seq);
	    gzclose(gzFP);
        if(DEBUG_COUNTING) fprintf(stderr,"Counting k-mers in FASTA file is complete. There are %ld valid lines\n",(long)validLines);
    }
    else if (manager->inputType == INPUT_LINES)
    {
	    while( fgets ( currentLine, MAX_CHARS_PER_LINE - 10,manager->inputFP) != NULL ) 
	    {
		    int linelen = strlen(currentLine);
		    if(streamOneString(manager->KWTree,&currentLine[0],linelen,&manager->substringCounts[0])!=EXIT_SUCCESS)
		    {
			    fclose(manager->inputFP);
			    return EXIT_FAILURE;
		    }	
	    }
    }
    else if(manager->inputType == INPUT_FILE)
    {
        fprintf(stderr,"To be done later\n"); //TBD
        return EXIT_FAILURE;
    }	
    else
    {
        fprintf(stderr,"Invalid input type: %d\n",manager->inputType);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;	
}

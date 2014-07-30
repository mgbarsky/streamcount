#include "count_kmers.h"

// Initialize KSEQ reader
KSEQ_INIT(gzFile, gzread)

int streamAndCountOneFile(KWTCounterManager *manager) {	
	//MT - set number of threads
    if(manager->numberOfThreads > 1)
        omp_set_num_threads(manager->numberOfThreads);

    fprintf (stderr, "Number of threads = %ld \n", (long) manager->numberOfThreads); 
	char currentLine[MAX_CHARS_PER_LINE];
	
    gzFile gzFP; 
    kseq_t* seq;
        
    SC_INT validLines=0; //counter of processed lines
    
    size_t PER_THREAD = 128; //128 lines are submitted to each thread job
    size_t N_BUFFER = MIN (manager->numberOfThreads,MAX_NUMBER_OF_THREADS)*PER_THREAD; //upto MAX_NUMBER_OF_THREADS threads max
    
    BufferCell *buffer;
    int buffer_index = 0;
    int result = EXIT_SUCCESS;

    //if multi-threaded - allocate buffer for a loop
    if(manager->numberOfThreads > 1) {
        buffer= malloc(N_BUFFER * sizeof(BufferCell));

        for(int i=0; i< N_BUFFER;i++) {
            buffer[i].seq = NULL;
            buffer[i].len = 0;
        }
    }   
    
	//now it depends on the type of the input file
    if(manager->inputType == INPUT_FASTA)  {        
        if(!( gzFP = gzdopen ( fileno(manager->inputFP) , "r" )))  {
		    fprintf(stderr,"Could not open input file as gz for reading\n");
		    return EXIT_FAILURE;
	    }

        seq = kseq_init(gzFP);
        
        if(manager->numberOfThreads == 1)  { //sequential - non-buffered processing
            while(kseq_read(seq) >= 0)   {  // read each sequence and process it
		        if ( streamOneString(manager->KWTree, seq->seq.s, seq->seq.l,&manager->substringCounts[0] ) != EXIT_SUCCESS )  
			        return EXIT_FAILURE;		        
                validLines++;
                if (validLines % 100000 == 0)
                   fprintf (stderr, "\rProcessed sequence %ld of the input file", (long) validLines );          
            }
        }
        else { //multi-threaded processing        
            // read sequences into a buffer 
            int done_reading = 0;
            while (!done_reading)  {
                if(kseq_read(seq) >= 0)   {
                    int newLen = seq->seq.l;
                    if (newLen > buffer[buffer_index].len) //re-allocate memory to hold the string                    
                        buffer[buffer_index].seq = malloc(newLen * sizeof(char));                    
                    buffer[buffer_index].len = seq->seq.l;
                    strcpy (buffer[buffer_index++].seq,seq->seq.s);
                    validLines++; 

                    if (validLines % 100000 == 0)
                        fprintf (stderr, "\rProcessed sequence %ld of the input file", (long) validLines );          
                }
                else
                    done_reading = 1;                   

                if (buffer_index == N_BUFFER || done_reading) {                
                    #pragma omp parallel for schedule(dynamic, PER_THREAD)                 
                    for(size_t i=0; i < N_BUFFER; ++i)       {
		                int thread_result = streamOneStringMT(manager->KWTree, buffer[i].seq, buffer[i].len,&manager->substringCounts[0] ); 
                        if (result == EXIT_SUCCESS)
                            result = thread_result;                                      
                    }  
                    
                    if (result != EXIT_SUCCESS) 
		                return EXIT_FAILURE;
	                   
                    buffer_index=0;
                }  
            }      
        }
       
        kseq_destroy(seq);
	    gzclose(gzFP);       
    }
    else if (manager->inputType == INPUT_LINES) {
            
        if(manager->numberOfThreads == 1) { //sequential - non-buffered processing
                 
            while( fgets ( currentLine, MAX_CHARS_PER_LINE - 10, manager->inputFP) != NULL ) {
		        int linelen = strlen(currentLine);
		        if(streamOneString(manager->KWTree,&currentLine[0],linelen,&manager->substringCounts[0])!=EXIT_SUCCESS)
			        return EXIT_FAILURE;
                validLines++;

                if (validLines % 100000 == 0)
                    fprintf (stderr, "\rProcessed sequence %ld of the input file", (long) validLines );                	        
	        }
        }
        else  {
            // read sequences into a buffer 
            int done_reading = 0;
            while (!done_reading)  {
                if(fgets ( currentLine, MAX_CHARS_PER_LINE - 10, manager->inputFP) != NULL )  {
                    int newLen = strlen(currentLine);
                    if (newLen > buffer[buffer_index].len) //re-allocate memory to hold the string                    
                        buffer[buffer_index].seq = malloc(newLen * sizeof(char));                    
                    buffer[buffer_index].len = newLen;
                    strcpy (buffer[buffer_index++].seq,currentLine);

                    validLines++;
                    if (validLines % 100000 == 0)
                        fprintf (stderr, "\rProcessed sequence %ld of the input file", (long) validLines );
                }
                else
                    done_reading = 1;                   

                if(buffer_index == N_BUFFER || done_reading) {                
                    #pragma omp parallel for schedule(dynamic, PER_THREAD)                 
                    for(size_t i=0; i < N_BUFFER; ++i)       {
		                int thread_result = streamOneStringMT(manager->KWTree, buffer[i].seq, buffer[i].len,&manager->substringCounts[0] ); 
                        if (result == EXIT_SUCCESS)
                            result = thread_result;                                      
                    }  
                    
                    if (result != EXIT_SUCCESS) 
		                return EXIT_FAILURE;
	                   
                    buffer_index=0;
                }                
                
                
            }  
        }
    }
    else if(manager->inputType == INPUT_FILE)   {
        fprintf(stderr,"To be done later. If your input is multi-line and each line is not extremely long (recommended shorter than 80 characters),\n");
        fprintf(stderr," you can convert it into FASTA format (by adding a description line) and process it in a FASTA mode. \n"); //TBD
        return EXIT_FAILURE;
    }	
    else  {
        fprintf(stderr,"Invalid input type: %d\n",manager->inputType);
        return EXIT_FAILURE;
    }
   
    if(manager->numberOfThreads > 1) {
        for(int i=0; i< N_BUFFER;i++)  //free buffer memory            
            free (buffer[i].seq);            
        free (buffer);
    } 

    if (validLines == 0) {
        fprintf ( stderr, "\rNo valid lines in the input file. Check that the specified input type is set correctly. \n" );
        return EXIT_FAILURE;
    }
    fprintf (stderr, "\rProcessed all %ld input sequences.          \n", (long)validLines );	
    		
    return EXIT_SUCCESS;	
}

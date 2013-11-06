/**
PREPARES - after this conversion - works on all systems, this only works for Unix with fts.h

This will convert input files into cleaner input files
renumber each file 0...N and record in FILE_NAMES
removes invalid lines
converts valid characters a,c,g,t to lower case
breaks larger files into 100MB chunks - can be redefined in MAX_FILE_SIZE_BYTES
**/

#include <sys/types.h>     
#include <sys/stat.h>
#include <err.h>
#include <fts.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#define MAX_CHARS_PER_LINE 10000
#define MAX_FILE_SIZE_BYTES 100000000 //to keep suffix arrays small and in memory - 100MB
#define MAX_NEW_FILE_NAME 100
#define MAX_PATH_LENGTH 250
#define DEBUG_MODE 0

static int  traverseFolder(char * const argv[]);

/*Input parameter: full path and name of a folder where input files are
Outputs: 
1. set of files inputs/input_ in folder inputs (which should be created before running), 
	each with valid DNA characters lines and of size no more than MAX_FILE_SIZE_BYTES
	files are numbered from 0 to numFiles-1
2. file inputs/FILE_NAMES where each line contains a real file name, which is mapped to a line number: 
	file previously names fileA becomes inputs/input_4 with fileA is written on the fifth line of file FILE_NAMES
*/
int
main(int argc, char * const argv[])
{
	int rc;

    	if ((rc = traverseFolder(argv + 1)) != 0)
        	rc = 1;
	return rc;
}

int prepareLine(char *from, char *to, int *size)
{
	int i=0;
	int endOfLine=0;
	while(!endOfLine)
	{
		if(from[i]==10 ||from[i]==32 ||from[i]=='\0')
		{
			endOfLine=1;
			*size=i;
			to[i]='\0';
			return 1;
		}
		else
		{
			char curr = from[i];
			
			switch (curr)
			{
				case 'a': case 'c': case 'g': case 't': case 'A': case 'C': case 'G': case 'T':					 				
					to[i]=toupper(curr);
					break;
				default:
					return 0; //invalid character encountered - ignore the line
			}
			i++;
		}		
	}
	return 1;
}

static int
traverseFolder(char * const argv[])
{
	FTS *ftsp;
	FTSENT *p, *chp;
	int fts_options = FTS_COMFOLLOW | FTS_LOGICAL | FTS_NOCHDIR;

	//*****************
	FILE *fpRead;
	FILE *fpWrite;
       
	FILE *fpFileNames;

	char current_line[MAX_CHARS_PER_LINE];
	char new_line [MAX_CHARS_PER_LINE];
	
	int currentFileID = 0;
	char currOutputName [MAX_PATH_LENGTH];
	int currentLineLength = 0;

	int lineResult = 0;

	long totalCharacters = 0;

	if ((ftsp = fts_open(argv, fts_options, NULL)) == NULL) {
		warn("fts_open");
		return -1;
	}

	// Initialize ftsp with as many argv[] parts as possible.
	chp = fts_children(ftsp, 0);
	if (chp == NULL) {
		return 0;               /* no files to traverse */
	}
	
	// Open file names file for writing - to map each input file name to its correponding file number
	if(!(fpFileNames= fopen ( "inputs/FILE_NAMES" , "w" )))
	{
		printf("Could not open FILE_NAMES file for writing \n");
		return (-1);
	}
	
	
	while ((p = fts_read(ftsp)) != NULL) 
	{
		if(p->fts_info == FTS_F) 
		{
			printf("processing file %d: %s\n", currentFileID, p->fts_name);   
			int filenamelength = strlen(p->fts_name);
			char lastchar=p->fts_name[filenamelength-1];
			if(DEBUG_MODE)
				printf("last char = %c\n",lastchar);
			
			// opening file for reading
			//ignore backup files
			if(lastchar!='~')
			{ 
				//open file for reading
				if(!(fpRead= fopen ( p->fts_path , "r" )))
				{
					printf("Could not open file \"%s\" for reading \n", p->fts_path);
					return (-1);
				}				

				//Record file name into INFO file
				fprintf(fpFileNames, "%s\n", p->fts_name);

				sprintf(currOutputName,"inputs/input_%d", currentFileID++);

				
				// opening file for writing
				if(!(fpWrite= fopen ( currOutputName , "w" )))
				{
					printf("Could not open file \"%s\" for writing \n", currOutputName);
					return (-1);
				}
				
				totalCharacters =0;
				//reading lines
				while( fgets (current_line, MAX_CHARS_PER_LINE-10, fpRead)!=NULL ) 
				{
					//this will replace valid characters a,c,g,t or A,C,G,T with lower case, 
					//and in case there are invalid characters will return false and this line will be ignored 
					lineResult = prepareLine(current_line,new_line,&currentLineLength);
				
					//writing concatenated content into a new file
					if(lineResult)
					{
						if(DEBUG_MODE)
							printf("%s\n",new_line);
						totalCharacters+=currentLineLength;
						if(totalCharacters<=MAX_FILE_SIZE_BYTES) //add to current file
							fprintf(fpWrite, "%s\n", new_line);
						else  //open a new file and write there
						{
							//close open output file
							fclose(fpWrite);
							//Record file name into INFO file - repeatingly
							fprintf(fpFileNames, "%s\n", p->fts_name);
							
							sprintf(currOutputName,"inputs/input_%d", currentFileID++);

				
							// opening new file for writing
							if(!(fpWrite= fopen ( currOutputName , "w" )))
							{
								printf("Could not open file \"%s\" for writing \n", currOutputName);
								return (-1);
							}
				
							totalCharacters =currentLineLength;
							fprintf(fpWrite, "%s\n", new_line);
						}		
					}				
				}
				fclose(fpRead);
				fclose(fpWrite);
			}
		}
	}
	fclose(fpFileNames);
	fts_close(ftsp);
	printf("Total %d new input files\n",currentFileID);	
	return 0;
}

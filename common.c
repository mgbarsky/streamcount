#include "common.h"

//printing data structures while in debug mode
void printPatterns (char **patterns, int numPatterns) {
	int i;
	fprintf(stderr,"\n k-mers:\n");
	for(i = 0; i < numPatterns; i++)	{
		fprintf(stderr,"%s\n",patterns[i]);
	}
}

//converts input line (while in a text mode) into a line of chars from a valid unified alphabet
//if the line is less than k chars it is considered invalid, if the mode is KMERS_FROM_LINES
int nextValidLineTextFile (FILE *inputFP, SC_INT k, char *outputLine) {
    char currentLine[MAX_CHARS_PER_LINE];
    if(fgets (currentLine, MAX_CHARS_PER_LINE-10, inputFP)!=NULL )   {
        int lineLen = lenValidChars(currentLine,strlen(currentLine));
		if(lineLen>=k) {            
            memcpy(outputLine,currentLine,lineLen);
            outputLine[lineLen]=0;
            return EXIT_SUCCESS;	    
        }
    }   
    return EXIT_FAILURE;
}

//combining counts for substrings obtained through kw tree and the mapping of each kmer (and its RC) to the corresponding kwtree leaf
//combinig these two data into counts for each k-mer
int combineSubstringCountsIntoKmersCounts(SC_INT  totalSubstringCounts,SC_INT *substringCounts,
    SC_INT totalKmers, KmerInfo *kmersInfo, 
    SC_INT *kmersCounts, int includeRC )
{
   
    SC_INT i, totalCount,countID,rcID;
    char repeated;

    //fprintf(stderr,"Total substring counted = %ld, total input kmers=%ld\n",(long)totalSubstringCounts,(long)totalKmers);
   
    //go in a loop through kmers info array and add total count for each k-mer - its count + rc count
	//if 'repeated' is true - set it to a negative value - because it will be ignored in subsequent calculations
	for(i=0; i< totalKmers; i++)
	{
		KmerInfo* kmer = &kmersInfo[i];
		countID = kmer->counterID;
		rcID = kmer->rcCounterID;
		repeated = kmer->repeated;
        
        if(includeRC)
        {
		    if(countID <= 0  || countID >= totalSubstringCounts || rcID==0 || rcID >= totalSubstringCounts)
		    {
				    fprintf(stderr,"Invalid value for %ld-th kmer-to-substring mapping: maps k-mer to %ld and rc to %ld\n",(long)i,(long)countID,(long)rcID);
				    return EXIT_FAILURE;
		    }

		    totalCount = substringCounts[countID] + substringCounts [rcID];
        }
        else
        {
            if(countID <= 0  || countID >= totalSubstringCounts )
		    {
				    fprintf(stderr,"Invalid value for %ld-th kmer-to-substring mapping: maps k-mer to %ld \n",(long)i,(long)countID);
				    return EXIT_FAILURE;
		    }

		    totalCount = substringCounts[countID];
        }
		if(repeated)
			totalCount=-totalCount;
		kmersCounts[i] = totalCount;
	}
    return EXIT_SUCCESS;
}

//Recursive printing of the tree.
//In order to print the keyword tree starting from the root call printKeywordTree(tree, 0,0)
void printKeywordTree (KWTNode *KWtree,SC_INT parentID, int level )
{
	SC_INT i,t;
	SC_INT newParentID;
	if(level == 0)
	{
		//print root
		fprintf(stderr,"Root\n");
	}
	
	//check if this is a leaf node
	if(KWtree[parentID].children[0] < 0) {
		//print indent
		for(t=0;t<level;t++) {
			fprintf(stderr," ");
		}		
		fprintf(stderr,"Leaf node <my id=%ld> for pattern %ld. Slink=%ld\n",(long)parentID, 
		                                        (long)(-KWtree[parentID].children[0]),
		                                        (long) KWtree[parentID].suffixLinkID);
	}
	
	for(i=0;i<SIGMA;i++) {
		if(KWtree[parentID].children[i] > 0) {
			//print indent
			for(t = 0; t < level; t++) {
				fprintf(stderr," ");
			}		
			fprintf(stderr,"Node %ld <parent id=%ld> <my id=%ld> Slink=%ld\n",(long)i,(long)parentID,
                    (long)KWtree[parentID].children[i],(long)KWtree[KWtree[parentID].children[i]].suffixLinkID);
			newParentID = KWtree[parentID].children[i];
			printKeywordTree (KWtree, newParentID, level+1);			
		}
	}
}

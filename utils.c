#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "KeywordTree.h"

int validChar(char c)
{
	switch( c ) 
	{
		case 'a': case 'A': 			
				return 1;
		case 'c': case 'C': 
				return 1;
		case 'g': case 'G': 
				return 1;
		case 't': case 'T': //3
				return 1;
		default :		
				return 0;
	}
}

int lenValidChars(char *currentLine, int lineLen)
{
	int i, end=0;
	int result=0;
	for(i=0;i<lineLen && !end;i++)
	{
		if(validChar(currentLine[i]))
			result++;
		else
			end=1;
	}
	return result;
}

INT getCharValue(char c)
{
	switch( c ) 
	{
		case 'a': case 'A': //0				
			return 0;
		case 'c': case 'C': //1
			return 1;
		case 'g': case 'G': //2
			return 2;
		case 't': case 'T': //3
			return 3;
		default :		
			printf("Invalid character encountered \n");
			exit(0);
	}
}

char getCharFromINT(INT n)
{
	switch( n ) 
	{
		case 0: //A				
				return 'A';
		case 1:  //C
				return 'C';
		case 2:  //G
				return 'G';
		case 3:  //T
				return 'T';
		default :		
				printf("Invalid number-to-character encountered \n");
				exit(0);
	}
}

INT getValidCharValue(char c)
{
	switch( c ) 
	{
		case 'a': case 'A': //0				
			return 0;
		case 'c': case 'C': //1
			return 1;
		case 'g': case 'G': //2
			return 2;
		case 't': case 'T': //3
			return 3;
		default :		
			return -1; //invalid char
	}
}

int validCharsToIntArray(char *line,int lineLength,INT *output)
{
	int i;
	int totalInts = 0;
	int end=0;
	for(i=0;i<lineLength && !end;i++)
	{
		INT curr = getValidCharValue(line[i]);
		if(curr==-1)
			end=1;
		else
		{
			output[totalInts++] = curr;
		}	
	}
	return totalInts;
}

//debug prints
void printPatterns(char **patterns, int numPatterns)
{
	int i;
	printf("\n UNIQUE PATTERNS:\n");
	for(i=0;i<numPatterns;i++)
	{
		printf("%s\n",patterns[i]);
	}
}

//to print tree starting from the root call printKeywordTree(tree, 0,0)
void printKeywordTree (KWTNode *KWtree,INT parentID, int level )
{
	int i,t;
	INT newParentID;
	if(level==0)
	{
		//print root
		printf("Root\n");
	}
	
	//check if this is a leaf node
	if(KWtree[parentID].children[0]<0)
	{
		//print indent
			for(t=0;t<level;t++)
			{
				printf(" ");
			}		
		printf("Leaf node <my id=%d> for pattern %d. Slink=%d\n",parentID, 
			-KWtree[parentID].children[0],
			KWtree[parentID].suffixLinkID);
	}
	
	for(i=0;i<SIGMA;i++)
	{
		if(KWtree[parentID].children[i]>0)
		{
			//print indent
			for(t=0;t<level;t++)
			{
				printf(" ");
			}		
			printf("Node %d <parent id=%d> <my id=%d> Slink=%d\n",i,parentID,KWtree[parentID].children[i],KWtree[KWtree[parentID].children[i]].suffixLinkID);
			newParentID = KWtree[parentID].children[i];
			printKeywordTree(KWtree,newParentID,level+1);			
		}
	}
}

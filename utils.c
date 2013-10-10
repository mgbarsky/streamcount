#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "KeywordTree.h"

int dictFromCharToInt[256]  = {
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1, 0,-1, 1,-1,-1,
-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1, 3,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1, 0,-1, 1,
-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1, 3,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1
};

char dictFromIntToChar[4]  ={'A','C','G','T'};
int validChar(char c)
{
	return dictFromCharToInt[(int)(c)]<0?0:1;
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

int getCharValue(char c)
{
	return dictFromCharToInt[(int)(c)];
}

char getCharFromINT(INT n)
{
	if(n<4)
		return 	dictFromIntToChar[n];
	else	
		printf("Invalid number-to-character encountered \n");
		exit(0);	
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

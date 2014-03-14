#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "General.h"

int dictFromCharToInt[256]  = {
-1,-1,-1,-1,-1,-1,-1,-1,-1,-2,
-2,-1,-2,-2,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,
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

char dictFromCharToComplement[256]  = {
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1, 'T',-1, 'G',-1,-1,
-1, 'C',-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1, 'A',-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1, 'T',-1, 'G',
-1,-1,-1, 'C',-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1, 'A',-1,-1,-1,
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

char getCharFromNumber(int n)
{
	if(n<4)
		return 	dictFromIntToChar[n];
	//else	
	fprintf(stderr,"Invalid number-to-character encountered \n");
	exit(0);	
}

char getComplement(char c)
{
	return dictFromCharToComplement[(int)c];
}

int produceReverseComplement(char *pattern, char *rcPattern)
{
	int i,r;
	int len=strlen(pattern);

	for(i=len-1, r=0; i>=0; i--,r++)
	{
		char complement=getComplement(pattern[i]);
		if(complement == -1)
		{
			fprintf(stderr,"INVALID CHARACTER ENCOUNTERED for complement: %s\n", pattern);
			return 1;
		}
		rcPattern[r] = complement;
	}
	rcPattern[r] ='\0';
	return 0;
}

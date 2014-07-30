#include "common.h"

/** Implementation of character convertions according to the DNA alphabet.
For speed, all conversions and validations are performed by accessing an entry in a constant array 
**/
const int dictFromCharToInt[256]  = {
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

const char dictFromCharToComplement[256]  = {
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

const char dictFromIntToChar[4]  ={'A','C','G','T'};

inline int validChar(char c){
	return dictFromCharToInt[(int)(c)] < 0 ? 0:1;
}

inline int lenValidChars(char *currentLine, int lineLen){
	int i, end=0;
	int result=0;
	for(i=0;i<lineLen && !end;i++)	{
		if(validChar(currentLine[i]))
			result++;
		else
			end=1;
	}
	return result;
}

inline int getCharValue(char c){
	return dictFromCharToInt[(int)(c)];
}

inline char getComplement(char c){
	return dictFromCharToComplement[(int)c];
}

inline int produceReverseComplement(char *pattern, char *rcPattern){
	int i,r;
	int len=strlen(pattern);

	for(i=len-1, r=0; i>=0; i--,r++){
		char complement=getComplement(pattern[i]);
		if(complement <0){
			fprintf(stderr,"INVALID CHARACTER ENCOUNTERED for complement: %s\n", pattern);
			return EXIT_FAILURE;
		}
		rcPattern[r] = complement;
	}
	rcPattern[r] ='\0';
	return EXIT_SUCCESS;
}

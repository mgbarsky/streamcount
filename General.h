#define MIN(a, b) ((a)<=(b) ? (a) : (b))
#define MAX(a,b) ((a)>=(b) ? (a) : (b))

#define MAX_PATH_LENGTH 500
#define MAX_CHARS_PER_LINE 100000  //make sure that patterns and input sequences are in lines no longer than that, or change it for longer lines 

#define snprintf snprintf 

//commonly used utils - dna-specific, implemented in 'dnautils.c'
int validChar(char c); //answers whether char c is from DNA alphabet
int lenValidChars(char *currentLine, int lineLen); //returns number of valid DNA chars in line

int getCharValue(char c); //returns 0,1,2,3 for a(A),c(C),g(G),t(T), returns -1 for all other chars
char getCharFromNumber(int n); //returns A,C,G,T for 0,1,2,3, terminates program if something else is requested
int produceReverseComplement(char *pattern, char *rcPattern); //produces char array rcPattern, complementary to char array pattern read in reverse order. Returns 1 if invalid (non-DNA character) encountered






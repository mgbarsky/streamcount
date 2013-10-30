#define MIN(a, b) ((a)<=(b) ? (a) : (b))
#define MAX(a, b) ((a)>=(b) ? (a) : (b))

#define INT int32_t  //this depends on the total length of the pattern set - which cannot exceed MAX(INT)
#define MAX_SET_SIZE 2000000000 //depends on the previous definition - maximum value of INT 
//final max size is minimum between this and number of nodes in the available memory
#define UINT uint32_t
#define MAX_COUNT 4000000000//depends on the previous definition of UINT - max value of UINT - to not to cause overflow
#define SIGMA 4  //alphabet size
#define MAX_PATH_LENGTH 500
#define MAX_CHARS_PER_LINE 100000  //make sure that patterns and input sequences are in lines no longer than that, or change it for longer lines 
#define MAX_K 100
#define DEBUG_KMERS 0
#define DEBUG_KWTREE 0
#define DEBUG_SERIALIZATION 0
#define DEBUG_COUNT 0


#define INPUT_LINES 0  //input for counting is provided as a file with lines. Max line length is MAX_CHARS_PER_LINE
#define INPUT_FILE 1 //input is a file (may contain multiple lines which should be considered as a single string). We need to consider the entire string running across multiple lines
#define INPUT_SNIPS 2 //input is provided as a context line around polymorphic site, and at the end of the line - comma-separated polymorphic nucleotides, which go directly into the middle of the line

typedef struct GlobalArgs {
    	char *patternFileName;      	/* -p option */
    	int memoryInMB;             	/* -m option */
    	int countOrNot;    		/* -c option */
    	int k;				/* -k option*/
    	int inputType;			/* -i option */
	int includeReverseComplement;	/* -r option*/	
	int isOutputDirectory;		/* -o option */
	char *outputDirName;
	int isInputDirectory;    	/* -d option */
    	char *inputDirName;          	/* directory to append before each input file */
    	int isFileWithFileNames;        /* -f option */
	char *fileFileNames;    	/* name of file containing file names */
	int numInputFiles;
	int inputFilesFromCmdLine;
	char **inputFiles; 		/* input files entered from the command line*/
} GlobalArgs;

//defines a node in the keyword tree
typedef struct KWTNode
{
	INT children[SIGMA]; //if 0 - no child, if negative at position 0 - leaf node - negated index in the array of counts
	INT suffixLinkID; //if negative - no link - only happens for the root node - all the rest have links	
}KWTNode;

//for breadth first traversal of the tree - to add suffix links
typedef struct Queue  
{
	int first;
	int freeSpot;
	int counter;
	INT *nodePointers;
}Queue;

//bookkeeping for using KWTree for search
typedef struct KWTCounterManager 
{
	INT treeSlotsNum; //number of slots in the array of tree nodes
	int k; //k - length of each k-mer
	INT totalPatterns; //total unique patterns to search - we are going to store that many counters per file, counters start from 1 in this array
	KWTNode *KWTree; //keyword tree holding all patterns to be counted
	UINT *patternCounts; //holds resulting count of each pattern - meaningful value begins from 1 in this array
	char inputFileName[MAX_PATH_LENGTH]; //the name of a current file where counting is performed
	char inputFileNamePrefix[MAX_PATH_LENGTH]; //deprecated left for compatibility with version 1
	int currentFileID; //deprecated left for compatibility with version 1	
}KWTCounterManager;

typedef struct KmerInfo 
{
	INT counterID; //index in an array of counters - starts from 1
	INT rcCounterID; //index in an array of counters of the reverse complement starts from 1		
	int startPosInLine; //starting position for this k-mer in the line of a pattern set file
	int lineNumber; //at what line does it start in a pattern set file
	char repeated; //how many times this k-mer occurs in the pattern set
}KmerInfo;

 //bookkeeping for building KWtree for pattern set - including repeating k-mers and reverse complement of each k-mer, which is counted as the same k-mer
typedef struct KWTreeBuildingManager
{
	INT maxSetSize;  //how much is allowed for a given memory
	INT totalSetSize;	
	INT treeSlotsNum; //number of slots in actual tree
	INT treeLeavesNum; //number of leaves in the tree eq. number of unique k-mers
	int k; //k - length of each k-mer
	INT totalPatterns; //number of total patterns, actual number of patterns is one less - since we start from 1
	INT maxPatterns; //to allocate memory
	int inputType;			//0-lines, 1 -file, 2-snips
	int includeReverseComplement;	//whether to include reverse complement	
	KWTNode *KWtree; //keyword tree holding all patterns to be counted
	KmerInfo *kmers; //string and mapping information about each k-mer
}KWTreeBuildingManager;

//information about previously saved KWtree
typedef struct KWTreeInfo 
{
	INT treeSlotsNum; //number of slots in an array of nodes
	int k; //k - length of each k-mer
	INT totalPatterns; //total unique patterns to search - we are going to store that many(-1) counters per file
}KWTreeInfo;

int process(GlobalArgs* globalArgs);
int buildPatternIndex(GlobalArgs *globalArgs);
int countAll(GlobalArgs *globalArgs);
//In file 'pattern_set_to_kwtree.c'
//fills in an array of k-mers from a given file
int fillPatternsArray(FILE *inputFP, char ** patterns, KmerInfo *patternsInfo, int64_t *totalPatterns, int k, int inputType);

//collects total number of non-unique k-mers - to allocate memory
int collectPatternsStats(FILE *inputFP, int k, int inputType, int64_t *totalPatterns); 

//preprocesses a set of patterns into a keyword tree of unique patterns
int preprocessPatternSet(KWTreeBuildingManager *manager, char *patternsFileName, int availableRamMB, int64_t* totalUniquePatterns); 

//in 'buildkwtree_main.c' 
int freeMemoryAfterKWtBuild (KWTreeBuildingManager* manager);

//in file 'keyword_tree.c'
int buildKeywordTree (KWTreeBuildingManager *manager, char **patterns, KmerInfo *patternsInfo, int64_t *numPatterns, int64_t *totalUniquePatterns);
int streamOneStringUnchanged(KWTCounterManager *manager, char *input, int strlength);
int traverseAndRecordPatterns(KWTreeBuildingManager *manager, char *currentPattern, int posInPattern,
	char **patterns, INT *newPatternCounter, int parentNodeID);
int addSuffixLinks (KWTNode *tree, int totalNodes);

//in file 'count_one_file.c' for streaming one input file - each line through the keyword tree
int streamAndCountOneFile(KWTCounterManager *manager);

//commonly used utils - in file 'utils.c'
int validChar(char c);
int lenValidChars(char *currentLine, int lineLen);
int numberKmersFromSnips(char *currentLine, int lineLen, int k) ; //parses the line to find out the number of polymorphic nucleotides
int getCharValue(char c);
char getCharFromINT(INT n);
int validCharsToIntArray(char *line,int lineLength,INT *output);
int produceReverseComplement(char *pattern, char *rcPattern);
int validSnip(char *line, int lineLen, int *snipPos, int *snipLen, int *charVal);

//debug prints
void printPatterns(char **patterns, int numPatterns);
void printKeywordTree (KWTNode *KWtree,INT parentID, int level );


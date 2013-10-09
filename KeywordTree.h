#define MIN(a, b) ((a)<=(b) ? (a) : (b))
#define MAX(a, b) ((a)>=(b) ? (a) : (b))

#define INT int32_t  //this depends on the total length of the pattern set - which cannot exceed MAX(INT)
#define MAX_SET_SIZE 2000000000 //depends on the previous definition - maximum value of INT 
//final max size is minimum between this and number of nodes in the available memory
#define UINT uint32_t
#define MAX_COUNT 4000000000//depends on the previous definition of UINT - max value of UINT - to not to cause overflow
#define SIGMA 4  //alphabet size
#define MAX_PATH_LENGTH 500
#define MAX_CHARS_PER_LINE 10000  //make sure that patterns and input sequences are in lines no longer than that, or change it for longer lines 

#define DEBUG_KWTREE 0
#define DEBUG_SERIALIZATION 0
#define DEBUG_COUNT 0

//defines a node in the keyword tree
typedef struct KWTNode
{
	INT children[SIGMA]; //if 0 - no child, if negative - leaf node - negated index in the array of counts
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
	char inputFileNamePrefix[MAX_PATH_LENGTH];
	int currentFileID;	
}KWTCounterManager;

 //bookkeeping for building KWtree for pattern set
typedef struct KWTreeBuildingManager
{
	INT maxSetSize;  //how much is allowed for a given memory
	INT totalSetSize;	
	KWTNode *KWtree; //keyword tree holding all patterns to be counted
	INT treeSlotsNum; //number of slots in actual tree
	int k; //k - length of each k-mer
	INT totalPatterns; //number of total patterns, actual number of patterns is one less 
	INT maxPatterns; //to allocate memory	
}KWTreeBuildingManager;

//information about previously saved KWtree
typedef struct KWTreeInfo 
{
	INT treeSlotsNum; //number of slots in an array of nodes
	int k; //k - length of each k-mer
	INT totalPatterns; //total unique patterns to search - we are going to store that many(-1) counters per file
}KWTreeInfo;

//In file 'pattern_set_to_kwtree.c'
//fills in an array of k-mers from a given file
int fillPatternsArray(FILE *inputFP, char **patterns, int64_t *totalPatterns, int k); 

//collects total number of non-unique k-mers - to allocate memory
int collectPatternsStats(FILE *inputFP, int k, int64_t *totalPatterns); 

//preprocesses a set of patterns into a keyword tree of unique patterns
int preprocessPatternSet(KWTreeBuildingManager *manager, char *patternsFileName, int availableRamMB, int64_t* totalUniquePatterns); 

//in 'buildkwtree_main.c' 
int freeMemoryAfterKWtBuild (KWTreeBuildingManager* manager);

//in file 'keyword_tree.c'
int buildKeywordTree (KWTreeBuildingManager *manager, char **patterns, int64_t *numPatterns, int64_t *totalUniquePatterns);
int streamOneString(KWTCounterManager *manager, INT *input, int length);
int traverseAndRecordPatterns(KWTreeBuildingManager *manager, char *currentPattern, int posInPattern,
	char **patterns, INT *newPatternCounter, int parentNodeID);
int addSuffixLinks (KWTNode *tree, int totalNodes);

//in file 'count_one_file.c' for streaming one input file - each line through the keyword tree
int streamAndCountOneFile(KWTCounterManager *manager);

//commonly used utils - in file 'utils.c'
int validChar(char c);
int lenValidChars(char *currentLine, int lineLen);
INT getCharValue(char c);
char getCharFromINT(INT n);
int validCharsToIntArray(char *line,int lineLength,INT *output);

//debug prints
void printPatterns(char **patterns, int numPatterns);
void printKeywordTree (KWTNode *KWtree,INT parentID, int level );


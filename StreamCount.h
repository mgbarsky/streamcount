#define SIGMA 4  //alphabet size
#define KMERS_FROM_LINES 0 //if extract k-mers from each line of the input
#define KMERS_FROM_FILE 1 //if extract k-mers from each position of the input file - spanning cross-lines

#define INPUT_FASTA 0
#define INPUT_FILE 2
#define INPUT_LINES 1

#define ALPHABET_TYPE_DNA 0

#define IGNORE_DUPLICATE_KMERS 0

#define INT int32_t //this depends on the total length of the pattern set - which cannot exceed MAX(INT). In case it does, replace this type with int64_t
#define MAX_COUNT 2000000000 //depends on the previous definition of INT - max value of INT - to not to cause overflow
#define MAX_INT 2000000000 //depends on the previous definition of INT - max value of INT - to not to cause overflow
//*****************
//Part 1. Preprocess set of k-mers into a keyword tree: read input according to its definition, extract k-mers and build KW tree with mapping
//********************

//defines a node in the keyword tree
typedef struct KWTNode
{
	INT children[SIGMA]; //if 0 - no child, if negative at position 0 - leaf node - negated index in the array of counts
	INT suffixLinkID; //if negative - no link - only happens for the root node - all the rest have links	
}KWTNode;

//for breadth first traversal of the tree - to add suffix links
typedef struct Queue  
{
	INT first;
	INT freeSpot;
	INT counter;
	INT *nodePointers;
}Queue;

typedef struct KmerInfo 
{
	INT counterID; //index in an array of counters - starts from 1
	INT rcCounterID; //index in an array of counters of the reverse complement starts from 1		
	INT startPosInLine; //starting position for this k-mer in the line of a pattern set file
	INT lineNumber; //at what line does it start in a pattern set file
	char repeated; //whether this k-mer is counted more than once
}KmerInfo;

//bookkeeping for building KWtree for pattern set - including repeating k-mers and reverse complement of each k-mer, which is counted as the same k-mer
typedef struct KWTreeBuildingManager
{
	long maxSetSize;  //how many tree nodes can be held in a given memory
	INT estimatedNumberOfLeaves; //estimated number of leaves - depends on include RC or not
	INT estimatedNumberOfKWTreeNodes;	//how many nodes it needs at most - to allocate memory
	INT actualNumberOfKWTreeNodes; //number of nodes in actual tree - known only after tree is built
	INT actualNumberOfLeaves; //corresponds to number of unique substrings contained in the tree
	INT maxNumberOfLeaves; //after reading actual patterns - know how many leaves can be in the tree
	INT treeLeavesNum; //number of leaves in the tree eq. number of UNIQUE k-mers
	INT k; //k - length of each k-mer
	INT originalNumberOfKmers; //number of total k-mers, actual number of patterns is one less - since we start from 1
    INT estimatedNumberOfKmers; //we need that in order to clean memory afterwards
	int inputType;			//0-lines, 1 -file
	int includeReverseComplement;	//whether to include reverse complement: default yes (1)
	KWTNode *KWtree; //keyword tree holding all patterns to be counted
	KmerInfo *kmersInfo; //string and mapping information about each k-mer
    char **kmers; //actual strings
}KWTreeBuildingManager;

//information about previously saved KWtree
/*typedef struct KWTreeInfo 
{
	INT numberOfTreeNodes; //number of slots in an array of nodes
	INT k; //k - length of each k-mer
	INT numberOfNonZeroLeaves; //total unique patterns to search - we are going to store that many(-1) counters per file
}KWTreeInfo;*/

//In file convertKmersIntoKWTree.c
int convertAllKmersIntoKWTreeReturnTree (FILE *kmersFP, int inputType, INT k, int includeRC, INT memoryMB, KWTreeBuildingManager *manager); //in-memory version
int fillKmersArrayAndInfo(KWTreeBuildingManager* manager, FILE *inputFP, 
	char **kmers, KmerInfo *kmersInfo, INT maxPossibleNumberOfKmers);
int collectKmerInputStats(FILE *inputFP, INT k, int inputType, INT *estimatedNumberOfKmers);
int convertAllKmersIntoKWTree (char *kmersFileName, int inputType, INT k, int includeRC, INT memoryMB);

//in file keyword_tree.c - all the operations on building and querying this data structure
int buildKeywordTree (KWTreeBuildingManager *manager, char ** patterns, KmerInfo *patternsInfo);
int addSuffixLinks (KWTNode *tree, int totalNodes);


//**************************
//Part 2 - Counting k-mers in one input file
//***************************
//bookkeeping for using KWTree for search
typedef struct KWTCounterManager 
{
	INT numberOfKWTreeNodes; //number of slots in the array of tree nodes
	INT k; //k - length of each k-mer
	INT numberOfKWTreeLeaves; //total unique patterns to search - we are going to store that many counters per file, counters start from 1 in this array
	KWTNode *KWTree; //keyword tree holding all patterns to be counted
	INT *substringCounts; //holds resulting count of each pattern - meaningful value begins from 1 in this array
	FILE *inputFP; //the pointer to a current file where counting is performed
    int inputType; //type of an input - default 0 - FASTA, LINES (1), FILE (2)
}KWTCounterManager;

int streamOneString(KWTNode* KWTree,char *input,int strlength,INT *patternCounts);
int streamAndCountOneFile(KWTCounterManager *manager);

//************************
//Part 3. Adding up counts for k-mer and its reverse complement
//************************
int combineSubstringCountsIntoKmersCounts(INT totalSubstringCounts,INT *substringCounts,
    INT totalKmers, KmerInfo *kmersInfo, 
    INT *kmersCounts, int includeRC );
int produceKmersCountsOneSample (char *countsBinaryFileName, char *mappingFileName, char *outputFileName, FILE *generalOutputFP, FILE *generalOutputFPText );
int nextValidLineTextFile (FILE *inputFP, INT k, char *outputLine);

//********************
//PRINTING DATA STRUCTURES WHILE IN DEBUG MODE
//*******************
#define DEBUG_KMERS_EXTRACTION 0
#define DEBUG_KWTREE 0
#define PRINT_KWTREE 0
#define DEBUG_COUNTING 0
void printPatterns(char **patterns, int numPatterns);
void printKeywordTree (KWTNode *KWtree,INT parentID, int level );


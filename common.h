#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>


//****************************
// TUNING STREAMCOUNT
// this values could be re-defined according to the application requirements
// 1. alphabet
// 2. number of k-mers to count
// 3. maximum number of threads - if you have more than 8 cores
//*****************************
#define SIGMA 4  //alphabet size - could be redefined for different alphabets

#define SC_INT int32_t //this depends on the total length of the pattern set - which cannot exceed MAX(SC_SC_INT). In case it does, replace this type with int64_t
#define MAX_COUNT 2000000000 //depends on the previous definition of SC_SC_INT - max value of SC_SC_INT - to not to cause overflow
#define MAX_SC_INT 2000000000 //depends on the previous definition of SC_SC_INT - max value of SC_SC_INT - to not to cause overflow

#define MAX_NUMBER_OF_THREADS 8
#define DEFAULT_NUMBER_OF_THREADS 6
//****************************
//DNA-specific functions and structs - can be re-defined to use streamcount for different alphabets
//*****************************

//produces char array rcPattern, complementary to char array *pattern, read in reverse order. Returns 1 if invalid (non-DNA character) encountered
int produceReverseComplement(char *pattern, char *rcPattern); 


typedef struct KmerInfo 
{
	SC_INT counterID; //index in an array of counters - starts from 1
	SC_INT rcCounterID; //index in an array of counters of the reverse complement starts from 1		
	SC_INT startPosInLine; //starting position for this k-mer in the line of a pattern set file
	SC_INT lineNumber; //at what line does it start in a pattern set file
	char repeated; //whether this k-mer occurs more than once in the entire k-mers file
}KmerInfo;
//********************
//END OF TUNING
//*******************


//type of input for k-mers
#define KMERS_FROM_LINES 0 //whether to extract k-mers from each line of the input
#define KMERS_FROM_FILE 1 //whether to extract k-mers from each position of the input file - spanning cross-lines

//type of input for files where to count k-mers
#define INPUT_FASTA 0 //k-mers are counted in each sequence of the FASTA or FASTQ file 
#define INPUT_LINES 1 //input file consists of text lines. k-mers are counted in each line.
#define INPUT_FILE 2 //input file consists of one large text, not divided into lines. k-mers are counted across lines in the entire file.
//Note: this input type is not currently supported. In order to simulate this behaviour, you can add a FASTA header to your file, and then 
//the entire text will be treated as one sequence.

//defines a node in the keyword tree used through all files
typedef struct KWTNode
{
	SC_INT children[SIGMA]; //if 0 - no child, if negative at position 0 - leaf node - negated index in the array of counts
	SC_INT suffixLinkID; //if negative - no link - only happens for the root node - all the rest have links	
}KWTNode;

//general macros
#define MIN(a, b) ((a)<=(b) ? (a) : (b))
#define MAX(a,b) ((a)>=(b) ? (a) : (b))

//general constants for allocating char arrays - strings
#define MAX_PATH_LENGTH 500
#define MAX_CHARS_PER_LINE 100000  //make sure that patterns and input sequences are in lines no longer than that, or change it for longer lines 

//this was for compatibility with windows version, which is not supported anymore
#ifdef _WIN32
	#define snprintf _snprintf
#endif

//********************
//PRINTING DATA STRUCTURES WHILE IN DEBUG MODE
//*******************
#define DEBUG_KMERS_EXTRACTION 0
#define DEBUG_KWTREE 0
#define PRINT_KWTREE 0
#define DEBUG_COUNTING 0
#define PRINT_COUNTING 0
void printPatterns(char **patterns, int numPatterns);
void printKeywordTree (KWTNode *KWtree,SC_INT parentID, int level );


//commonly used alphabet evaluation utils - declarations are not dna-specific. DNA-specific implementations are in 'dna_common.c'
//the idea is that you need to include common.h with any alphabet, but implement these declared functions in the alphabet-specific file
int validChar(char c); //answers whether char c is from a given alphabet
int lenValidChars(char *currentLine, int lineLen); //returns number of consecutive valid alphabet chars in line

int getCharValue(char c); //for a valid aphabet: for example for DNA returns 0,1,2,3 for a(A),c(C),g(G),t(T), returns -1 for all other chars

//converts input line (while in a text mode) into a line of chars from a valid unified alphabet
int nextValidLineTextFile (FILE *inputFP, SC_INT k, char *outputLine);

//this is designed to allow to write an exit status of each run of streamcount (success or failure) into a file
//so when performing parallel processing on cluster you know when all jobs are terminated and whether they were successful
int endProgram(int exitStatus, int indicateEnd, FILE *endFile );

#endif


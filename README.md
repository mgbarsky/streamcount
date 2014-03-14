streamcount
===========
This is a program which counts occurences of k-mers (strings of length k characters) in an arbitrarily large input.

The program first takes a set of pattern strings, breaks the strings into k-mers, and builds from this set of k-mers
a keyword tree with suffix links (see Aho-Corasick algorithm). 

In the second part, each line of an input file is streamed through the keyword tree 
and the counters of the corresponding k-mers for this file are collected.

The number of k-mers which can be simultaneously counted is limited by the amount of the available RAM.
The number is also limited by the use of a signed integer INT defined as int32_t on lines 13, 15 of StreamCount.h.
With this definition, we can build an index for at most Int32.MaxValue/k input k-mers.
To increase this limit, redefine INT as int64_t and recompile.

*************
Dependencies:
*************
zlib

*************
To compile:
*************
make

************
To run:
************
If you add a path to the compiled streamcount to your PATH variable, 
it can be run as a standard unix command: streamcount

To run a program (./streamcount) you need to specify the following parameters
****************
PROGRAM ARGUMENTS
****************

Mandatory:
**********

--kmers 'kmers_file'

where 'kmers_file' is the full path and file name of the file from which to extract the k-mers.
The file with k-mers should contain only characters from a valid DNA alphabet. 
This should be dealt with prior to running the program.

-i --input 'input_file'

where 'input_file' is the full path and file name of the file where to count the k-mers.


If the input option is not specified, the program tries to read the input text from stdin.
In this case, the following commands are valid:

cat 'input_file' |./streamcount --kmers 'kmers_file'
./streamcount --kmers 'kmers_file' < 'input_file'

By specifying only these two mandatory parameters, we accept the following default program behaviour:
1. 'input_file' is of type FASTA. It can be compressed.

2. Each line of 'kmers_file' is treated as a separate k-mer.

3. The final count for each k-mer includes a count for its reverse complement string.

4. The final counts for each k-mer are written to stdout, one count per line.

5. If the k-mers in 'kmers_file' are not unique, the information about this is supressed.

Optional:
*********
To modify default behavior:

Input options: 
**************
-k='k' 
length of each k-mer. If there are more than one k-mer in each input line, all of them will be considered. In this case, output for each line will consist of a line of comma-separated counts

--kmers-multiline
extract k-mers from 'kmers_file' treating the entire file as one string

--input-plain-text
treat input as text lines, rather than FASTA.
 
Counting options:
***************** 
--no-rc 
do not include count of reverse complement into final count of each k-mer. This option can be useful when counting k-mers in a genomic sequence.

-m,     --mem='MEMORY_MB'
amount of memory available for indexing k-mers. Specify the amount of memory (in MB) that you are ready to sacrifice to hold a k-mer index. 
This is used to estimate if you can hold k-mers index prior to processing. 
Default: 4000MB

Output options: 
*************** 
--printseq
print each original line of 'kmers_file' before its count(s). 

--repeat-mask-tofile='repeat-mask-file'
for each k-mer prints to 'repeat-mask-file' 0 or 1. 
1 is printed if this k-mer is not unique (repeats) in the 'kmers_file'
    



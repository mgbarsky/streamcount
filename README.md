streamcount
===========
This is a program which counts occurences of k-mers (strings of length k characters) in an arbitrarily large input.

The program first takes a set of pattern strings, breaks them into k-mers, and builds from this set of patterns 
a keyword tree with suffix links (see Aho-Corasick algorithm). This keyword tree is serialized to disk, 
to stream any input file against it.

In the second part, each line of each input file is streamed through the keyword tree 
and the counters of the corresponding k-mers for this file are collected and serialized to disk.

*************
To compile:
*************
make

************
To run:
************
To run a program (./streamcount) you need to specify the following parameters


Mandatory:
**********

-p 'file name'

where 'file name' is the full path and file name of the file from which to extract the k-mers

-k value of k

If only these parameters are specified, the program will create a keyword tree to be used later.
The output keyword tree will be written into the same directory as the 'file name' in form of 3 files:

1. 'file name'_k-mers_MAPPINGS - this is a binary file which contains an information for each extracted k-mer: 
line number, start position in this line, the position of the counter for this k-mer in the array of counters, 
the position of the counter for the reverse complement of this k-mer in the array of counters, and 0/1 flag 
which specifies whether this k-mer was unique in the input patterns file (0), or it occurs in more than 1 place

2. 'file name'_k-mers_KWTREE_INFO - binary file containing parameters of the keyword tree - to be loaded later

3. 'file name'_k-mers_KWTREE - binary file containing the tree itself

Note: if the files with the same name already exist, the tree will not be re-created.


Optional:
*********
-r If this option is specified, the reverse complement of each k-mer will be added to the set of k-mers

-i Input mode for the extraction of k-mers: 

  0 - each line in the pattern file is treated as a separate line;
  1 - all lines of the input file are treated as one long string

-c If this option is specified, program performs the k-mer counting in a set of input files.
In this case, you have to specify the input files where to perform the counting.

There are 2 options:

1. List all input files after all options, separated by space

2. Provide a file with a list of file names, where each line contains the file name to be processed: -f 'file name'

Additionaly, for both cases you can specify the relative path to the input files using option -d 'directory name'

Counting output:
The counting will be performed in each specified file, and an array with counters will be saved to disk in a binary file
'input fle name'_k-mers_COUNTS, in the same directory where the input files are.

If you want to write the counts into a different directory, you will need to specify an additional option
-o 'output directory'

Additional utils
****************
In this version there are two additional programs which convert binary outputs into text format, 
so that they can be post-processed using different file manipulation tools

./countstotext 'patterns file name fullpath' 'size of k-mers' 'input file name(s) full path'

This will convert counters array into a text file, each counter on a separate line

./patternsstotext 'patterns file name fullpath' 'size of k-mers'

This will convert the information about the k-mers into a tabular format with the following columns:
line number, start pos in this line, index in the array of counters, same index for reverse complement, repeated (no, yes)


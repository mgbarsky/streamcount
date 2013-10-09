streamcount
===========
This is a program which counts occurences of k-mers (aka strings of length k characters) in an arbitrarily large input.

The program first takes a set of pattern strings, breaks them into k-mers, and builds from this set of patterns 
a keyword tree with suffix links (see Aho-Corasick algorithm). This keyword tree is serialized to disk, 
to stream any input file against it.

In the second part, each line of each input file is streamed through the keyword tree 
and the counters of the corresponding k-mers for this file are collected and serialized to disk.


To compile:
make

To run:
*****************
Part 1
*****************
./buildkwtree 'full path and file name of the file with pattern set' 'value of k' 'memory available in MB'

Output: 
    binary file with a keyword tree of all k-mers
    a binary info file holding information about the tree
    MAPPING file with all unique k-mers found in the input: each unique k-mer is mapped to its corresponding line number
            meaningful k-mers start from line 1 (zero-line is ignored)
    All this files are in the same directory where the pattern file is

*****************
Part 2
*****************
If you need to count k-mers in a single file run:
./streamandcountonefile 'full path and name of the patterns file'  'full path and prefix of numbered input files' 
          'input file index' 'value of k'
The program assumes that you have preprocessed your input into a set of files with a shared prefix and a numeric index.
The preprocessing routine for a single directory (for Linux) can be found here: 
  https://github.com/mgbarsky/digest/blob/master/prepare_inputs.c

Output: for each unique k-mer in a MAPPING file - an array of counters, serialized in a binary file.
        Such an array is created and written separately for each input file
        The count files are in the same directory where the input files are

<h1>streamcount</h1>
This is a program which counts occurences of k-mers (strings of length k characters) 
in an arbitrarily large input.

The program first takes a set of pattern strings, breaks the strings into k-mers, 
and builds from this set of k-mers a keyword tree with suffix links (see Aho-Corasick algorithm). 

In the second part, each line of an input file is streamed through the keyword tree 
and the counters of the corresponding k-mers for this file are collected.

The number of k-mers which can be simultaneously counted is limited by the amount of the available RAM.
The number is also limited by the use of the signed integer SC_INT defined as int32_t on line 18 of common.h.
With this definition, we can build an index for at most Int32.MaxValue/k input k-mers.
To increase this limit, redefine SC_INT as int64_t and recompile.


<h2>Dependencies:</h2>
<pre> <code>zlib</code> </pre>

<h2>To compile:</h2>
<pre> <code>make</code> </pre>

<h2>To run:</h2>
th to the compiled streamcount to your PATH variable, 
it can be run as a standard unix command: streamcount

<h2>Program arguments</h2>

<h3>Required:</h3>
<pre> <code>--kmers 'kmers_file'</code> </pre>
where 'kmers_file' is the full path and file name of the file from which to extract the k-mers.
<br>NOTE: The file with k-mers should contain only characters from a valid DNA alphabet. 
This should be dealt with prior to running the program.

<pre> <code>-i --input 'input_file'</code> </pre>
where 'input_file' is the full path and file name of the file where to count the k-mers.
If the input option is not specified, the program tries to read the input text from stdin.
In this case, the following commands are valid:

<pre> <code>cat 'input_file' |./streamcount --kmers 'kmers_file'</code> </pre>
<pre> <code>./streamcount --kmers 'kmers_file' < 'input_file'</code> </pre>

By specifying only these two parameters, we accept the following default program behaviour:
<ol>
<li>'input_file' is of type FASTA. It can be compressed.</li>
<li>Each line of 'kmers_file' is treated as a separate k-mer.</li>
<li>The final count for each k-mer includes a count for its reverse complement string.</li>
<li>The final counts for each k-mer are written to stdout, one count per line.</li>
<li>If some k-mers in 'kmers_file' are not unique, the information about this is supressed.</li>
<li>Multi-threaded execution with DEFAULT_NUMBER_OF_THREADS defined on line 24 in common.h.</li>
</ol>

<h3>Optional:</h3>

<h4>Input options:</h4>
<pre> <code>-k='k'</code> </pre>
length of each k-mer. 
If there are more than one k-mer in each input line, all of them will be considered. 
In this case, output for each line will consist of a line of comma-separated counts

<pre> <code>--kmers-multiline</code> </pre>
extract k-mers from 'kmers_file' treating the entire file as one string

<pre> <code>--input-plain-text</code> </pre>
treat input as text lines, rather than FASTA.

<pre> <code>--t</code> </pre>
number of threads for multi-threaded processing. 
It is optimal to define the number of threads as the number of cores. 
Maximum number of threads is set to 8. It can be redefined in common.h line 23 
 
<h4>Counting options:</h4>
<pre> <code>--no-rc</code> </pre> 
do not include count of reverse complement into final count of each k-mer. 
This option can be useful when counting k-mers in a genomic sequence.

<pre> <code>-m,     --mem='MEMORY_MB'</code> </pre>
amount of memory available for indexing k-mers. 
Specify the amount of memory (in MB) that you are ready to sacrifice to hold a k-mer index. 
This is used to estimate if you can hold k-mers index prior to processing. 
Default: 4000MB

<h4>Output options:</h4>
<pre> <code>--printseq</code> </pre>
print each original line of 'kmers_file' before its count(s). 

<pre> <code>--repeat-mask-tofile='repeat-mask-file'</code> </pre>
for each k-mer prints to 'repeat-mask-file' 0 or 1. 
1 is printed if this k-mer is not unique (repeats) in the 'kmers_file'.
This is used if you need a precise count for all k-mers extracted from the same line. 
Because the same k-mer occurs also on a different line, the counts of consecutive k-mers are distorted.

<h2>Sample usage:</h2>
In folder sample_data.zip there are one sample input file, and one k-mers file.
Folder also contains SAMPLE_RUNS.txt with examples of running streamcount.



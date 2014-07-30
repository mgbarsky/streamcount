CC = gcc
CFLAGOPT = -O3 -Wall 
CFLAGS = -D_LARGEFILE_SOURCE
CFLAGS += -fno-exceptions
CFLAGS += -finline-functions
CFLAGS += -funroll-loops
CFLAGOFFSET = -D_FILE_OFFSET_BITS=64
LDFLAGS=-lz
MATHFLAG=-lm
MTFLAG=-fopenmp
MTFLAG += -std=c99

# Source files
SC_SRC=common.c dna_common.c keyword_tree.c kmers_to_kwtree.c count_kmers.c streamcount.c

# Binaries
all: streamcount

#streams the lines of the input file (in any format - fasta, text, compressed) and counts k-mers
streamcount: $(SC_SRC)
	$(CC) $(CFLAGOPT) $(CFLAGOFFSET) $(CFLAGS) $(MTFLAG) $^ -o $@ $(LDFLAGS)

clean:  
	rm streamcount

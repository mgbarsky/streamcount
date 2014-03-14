CC = gcc
CFLAGOPT = -O3 -Wall 
CFLAGS = -D_LARGEFILE_SOURCE
CFLAGS += -fno-exceptions
CFLAGS += -finline-functions
CFLAGS += -funroll-loops
CFLAGOFFSET = -D_FILE_OFFSET_BITS=64
LDFLAGS=-lz
MATHFLAG=-lm

# Source files
SC_SRC=generalutils.c dnautils.c keyword_tree.c convertKmersIntoKWTree.c countKmersInFile.c MAIN.c

# Binaries
all: streamcount 

#streams the lines of the input file (in any format - fasta, text, compressed) and counts k-mers
streamcount: $(SC_SRC)
	$(CC) $(CFLAGOPT) $(CFLAGOFFSET) $(CFLAGS) $^ -o $@ $(LDFLAGS)
clean:  
	rm streamcount 

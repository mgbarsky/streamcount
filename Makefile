CC = gcc
CFLAGOPT = -O3 -Wall 
CFLAGS = -D_LARGEFILE_SOURCE
CFLAGS += -fno-exceptions
CFLAGS += -finline-functions
CFLAGS += -funroll-loops
CFLAGOFFSET = -D_FILE_OFFSET_BITS=64
LDFLAGS=-lz

# Source files
SC_SRC=utils.c keyword_tree.c pattern_set_to_kwtree.c streamandcount_all.c main.c
CT_SRC=counters_binary_to_text.c
PT_SRC=kmers_binary_to_text.c

# Binaries
all: streamcount countstotext patternstotext

streamcount: $(SC_SRC)
	$(CC) $(CFLAGOPT) $(CFLAGOFFSET) $(CFLAGS) $^ -o $@ $(LDFLAGS)
countstotext: $(CT_SRC)
	$(CC) $(CFLAGOPT) $(CFLAGOFFSET) $(CFLAGS) $^ -o $@ $(LDFLAGS)
patternstotext: $(PT_SRC)
	$(CC) $(CFLAGOPT) $(CFLAGOFFSET) $(CFLAGS) $^ -o $@ $(LDFLAGS)
clean:  
	rm streamcount countstotext patternstotext

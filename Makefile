CC = gcc
CFLAGOPT = -O3 -Wall -pg
CFLAGS = -D_LARGEFILE_SOURCE
CFLAGS += -fno-exceptions
CFLAGS += -finline-functions
CFLAGS += -funroll-loops
CFLAGOFFSET = -D_FILE_OFFSET_BITS=64

all: streamcount countstotext patternstotext
streamcount:
	$(CC) $(CFLAGOPT) $(CFLAGOFFSET) $(CFLAGS)  utils.c keyword_tree.c pattern_set_to_kwtree.c streamandcount_all.c main.c -o streamcount
countstotext:
	$(CC) $(CFLAGOPT) $(CFLAGOFFSET) $(CFLAGS)  counters_binary_to_text.c -o countstotext
patternstotext:
	$(CC) $(CFLAGOPT) $(CFLAGOFFSET) $(CFLAGS)  kmers_binary_to_text.c -o patternstotext
clean:  
	rm streamcount countstotext patternstotext

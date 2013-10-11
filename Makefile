CC = gcc
CFLAGOPT = -O3 -Wall 
CFLAGS = -D_LARGEFILE_SOURCE
CFLAGS += -fno-exceptions
CFLAGS += -finline-functions
CFLAGS += -funroll-loops
CFLAGOFFSET = -D_FILE_OFFSET_BITS=64

all: streamcount
buildkwtree:
	$(CC) $(CFLAGOPT) $(CFLAGOFFSET) $(CFLAGS) utils.c pattern_set_to_kwtree.c keyword_tree.c buildkwtree_main.c -o buildkwtree
streamandcountonefile:
	$(CC) $(CFLAGOPT) $(CFLAGOFFSET) $(CFLAGS) utils.c keyword_tree.c count_one_file.c streamandcount_onefile_main.c -o streamandcountonefile
streamandcountall:
	$(CC) $(CFLAGOPT) $(CFLAGOFFSET) $(CFLAGS)  utils.c keyword_tree.c count_one_file.c streamandcount_all_main.c -o streamandcountall
streamcount:
	$(CC) $(CFLAGOPT) $(CFLAGOFFSET) $(CFLAGS)  utils.c keyword_tree.c pattern_set_to_kwtree.c streamandcount_all.c buildkwtree.c main.c -o streamcount
prepareonunix:
	$(CC) $(CFLAGOPT) $(CFLAGS)  prepare_inputs.c -o prepare
clean:  
	rm streamcount

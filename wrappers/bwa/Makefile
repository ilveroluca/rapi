

CC = gcc

WRAP_MALLOC = -DUSE_MALLOC_WRAPPERS
CFLAGS = -g -Wall -O2 -Wno-unused-function -std=c99
DFLAGS = -DHAVE_PTHREAD $(WRAP_MALLOC)
LIBS = -lm -lz -lpthread

# these depend on BWA_PATH
INCLUDES := -I../../ -I$(BWA_PATH)
SOURCES := rapi_bwa.c rapi_test.c $(BWA_PATH)/bwamem.c
# we define the object names from the source names, substituting the
# extension and keeping online the file name (removing the directory part)
OBJS := $(notdir $(SOURCES:.c=.o))

.SUFFIXES:.c .o

.PHONY: bwa clean

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDES) $(DFLAGS) $< -o $@

all: rapi

rapi: bwa $(BWA_PATH)/libbwa.a $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o rapi -L$(BWA_PATH) -lbwa $(LIBS) 

bwa:
	@echo "BWA_PATH is $(BWA_PATH)"
	$(if $(BWA_PATH),, $(error "You need to set the BWA_PATH variable on the cmd line to point to the compiled BWA source code (e.g., make BWA_PATH=/tmp/bwa)"))

bwamem.o: bwa $(BWA_PATH)/bwamem.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $(DFLAGS) $(BWA_PATH)/bwamem.c -o $@

clean:
	rm -f $(OBJS)

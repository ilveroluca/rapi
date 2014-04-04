
CC = gcc

WRAP_MALLOC = -DUSE_MALLOC_WRAPPERS
CFLAGS = -g -Wall -O2 -Wno-unused-function -std=c99
DFLAGS = -DHAVE_PTHREAD $(WRAP_MALLOC)
LIBS = -lm -lz -lpthread
INCLUDES= aligner.h

SOURCES = rapi_bwa.c rapi.c bwamem.c
OBJS = $(SOURCES:.c=.o)

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $< -o $@

all: rapi

rapi: libbwa.a $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o rapi -L. -lbwa $(LIBS) 

clean:
	rm -f $(OBJS)

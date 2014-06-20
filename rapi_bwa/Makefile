
CC = gcc
AR = ar

WRAP_MALLOC := -DUSE_MALLOC_WRAPPERS
CFLAGS := -g -Wall -O2 -Wno-unused-function -std=c99 -fPIC
DFLAGS := -DHAVE_PTHREAD $(WRAP_MALLOC)
LIBS := -lm -lz -lpthread
RAPI_LIB := librapi_bwa.a

# the includes depend on BWA_PATH
INCLUDES := -I../include/ -I$(BWA_PATH)

SOURCES := $(wildcard *.c)
# we define the object names from the source names, substituting the
# extension and keeping online the file name (removing the directory part)
OBJS := $(notdir $(SOURCES:.c=.o))

.SUFFIXES:.c .o

.PHONY: bwa clean

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDES) $(DFLAGS) $< -o $@

all: $(RAPI_LIB)

$(RAPI_LIB): $(OBJS)
	$(AR) -csru $@ $(OBJS)

bwa:
	@echo "BWA_PATH is $(BWA_PATH)"
	$(if $(BWA_PATH),, $(error "You need to set the BWA_PATH variable on the cmd line to point to the BWA source code (e.g., make BWA_PATH=/tmp/bwa)"))

clean:
	rm -f $(OBJS) $(RAPI_LIB)

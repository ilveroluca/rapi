
CC := gcc

WRAP_MALLOC := -DUSE_MALLOC_WRAPPERS
CFLAGS := -g -Wall -std=c99 -O2
DFLAGS := -DHAVE_PTHREAD $(WRAP_MALLOC)
LIBS := -lm -lz -lpthread

# the includes depend on BWA_PATH
INCLUDES := -I../include/

SOURCES := rapi_test.c
# we define the object names from the source names, substituting the
# extension and keeping online the file name (removing the directory part)
OBJS := $(notdir $(SOURCES:.c=.o))
EXE := rapi_test
RAPI_LIB := ../rapi_bwa/librapi_bwa.a

.SUFFIXES:.c .o

.PHONY: clean

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDES) $(DFLAGS) $< -o $@

all: $(EXE)

$(EXE): bwa $(BWA_PATH)/libbwa.a $(RAPI_LIB) $(OBJS)
	# we need to link to our static librapi_bwa.a as well as BWA's static library
	$(CC) $(CFLAGS) $(OBJS) -o $(EXE) -L$(BWA_PATH) -L$(dir $(RAPI_LIB)) -lrapi_bwa -lbwa $(LIBS) 

bwa:
	@echo "BWA_PATH is $(BWA_PATH)"
	$(if $(BWA_PATH),, $(error "You need to set the BWA_PATH variable on the cmd line to point to the compiled BWA source code (e.g., make BWA_PATH=/tmp/bwa)"))


clean:
	rm -f $(OBJS) $(EXE)

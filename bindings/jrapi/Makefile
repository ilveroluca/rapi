

Package := it.crs4.rapi.lowrapi

CC = gcc

SHELL=/bin/bash
Platform := $(shell uname | tr A-Z a-z)

JAVA_INCLUDES := -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/$(Platform)

# for dependencies
JavaFiles = $(shell find $(GENSRC) $(JAVASRC) -type f)

INTERFACE := jrapi.i
ROOT := ../..
INCLUDES := -I$(ROOT)/include $(JAVA_INCLUDES)
CFLAGS := -g -Wall -std=c99 -fPIC
# with optimizations, turn off strict-aliasing as per http://www.swig.org/Doc3.0/Java.html
CFLAGS := $(CFLAGS) -O2 -fno-strict-aliasing

RAPI_BWA := $(ROOT)/rapi_bwa/librapi_bwa.a
SHARED := jrapi.so
JAVASRC := src
GENSRC := gensrc
JAR := build/jrapi.jar


# I've seen swig installed under various names.  Look for all of them
SWIG := $(shell type -P swig swig3 swig3.{0..9} | tail  -n 1)

.SUFFIXES:.c .o
.PHONY: clean

all: jar

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDES) $(DFLAGS) $< -o $@

jrapi_wrap.c $(GENSRC): $(ROOT)/include/rapi.h $(INTERFACE)
	$(info "Using SWIG executable $(SWIG)")
	mkdir -p $(GENSRC)
	$(SWIG) -package $(Package) -java $(INCLUDES) -outdir $(GENSRC) $(INTERFACE)

#jrapi_wrap.o: jrapi_wrap.c
#	$(CC) $(CFLAGS) $(INCLUDES) -c jrapi_wrap.c

$(SHARED): jrapi_wrap.o
# to build the shared library we link against the static rapi_lib
ifeq ($(Platform), darwin)
	# LP:  no idea whether this works.  If you have a Mac test it out.
	$(CC) -shared $(CFLAGS) -o $@ jrapi_wrap.o -L$(dir $(RAPI_BWA)) -lrapi_bwa -lz
	install_name_tool -change $(SHARED)
else
	$(CC) -shared $(CFLAGS) -o $@ $< -L$(dir $(RAPI_BWA)) -lrapi_bwa -lz
endif

jar: $(JAR)

$(JAR): $(SHARED) $(JavaFiles)
	@echo "JavaFiles: $(JavaFiles)"
	ant package

clean:
	rm -f $(SHARED) jrapi_wrap.[co]
	rm -rf $(GENSRC)
	ant clean

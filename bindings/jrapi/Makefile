###################################################################################### 
#  This code is published under the The MIT License.
# 
#  Copyright (c) 2016 Center for Advanced Studies,
#                       Research and Development in Sardinia (CRS4), Pula, Italy.
# 
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
# 
#  The above copyright notice and this permission notice shall be included in all
#  copies or substantial portions of the Software.
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#  SOFTWARE.
###################################################################################### 



Package := it.crs4.rapi

CC = gcc

SHELL=/bin/bash
Platform := $(shell uname | tr A-Z a-z)

JAVA_INCLUDES := -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/$(Platform)

# for dependencies
JavaFiles = $(shell find $(GENSRC) $(JAVASRC) -type f -name '*.java')

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

debug: CFLAGS += -g -O0
debug: jar

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

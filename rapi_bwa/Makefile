
###############################################################################
# Copyright (c) 2014-2016 Center for Advanced Studies,
#                         Research and Development in Sardinia (CRS4)
# 
# Licensed under the terms of the MIT License (see LICENSE file included with the
# project).
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
###############################################################################

CC = gcc
AR = ar

WRAP_MALLOC := -DUSE_MALLOC_WRAPPERS
CFLAGS := -g -Wall -O2 -std=c99 -fPIC
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

bwa_header.h: bwa
	./extract_bwa_header.py $(BWA_PATH) > bwa_header.h

# We create the static rapi_bwa library ("RAPI_LIB") as a thin ar archive that
# also indexes the symbols in the static libbwa.a.  In this way we can link client
# only to our static library, without having to keep track of where libbwa.a is.
#
# XXX: Note that our *thin* archive doesn't contain the bwa code but merely
# references the static archive, so when the client code is built libbwa.a must
# still be in its original location.
$(RAPI_LIB): $(OBJS) $(BWA_PATH)/libbwa.a
	rm -f $@
	$(AR) -crT $@ $(OBJS) $(BWA_PATH)/libbwa.a

$(OBJS): bwa bwa_header.h

bwa:
	@echo "BWA_PATH is $(BWA_PATH)"
	$(if $(BWA_PATH),, $(error "You need to set the BWA_PATH variable on the cmd line to point to the BWA source code (e.g., make BWA_PATH=/tmp/bwa)"))

clean:
	rm -f $(OBJS) $(RAPI_LIB)

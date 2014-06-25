
CC = gcc

PYTHON_INCLUDES := $(shell python-config --includes)
INCLUDES := -I../include $(PYTHON_INCLUDES)
CFLAGS := -g -Wall -O2 -Wno-unused-function -std=c99 -fPIC
RAPI_BWA := ../rapi_bwa/librapi_bwa.a
SHARED := _rapi.so

.SUFFIXES:.c .o
.PHONY: bwa clean

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDES) $(DFLAGS) $< -o $@

all: $(SHARED)

rapi_wrap.c: ../include/rapi.h ../include/rapi.i
	swig -builtin -python $(INCLUDES) -o rapi_wrap.c ../include/rapi.i

rapi_wrap.o: rapi_wrap.c
	$(CC) $(CFLAGS) $(INCLUDES) -c rapi_wrap.c

$(SHARED): bwa rapi_wrap.o
	# to build the shared library, we link against our static rapi_lib and the static BWA library, libbwa.a
	$(CC) -shared $(CFLAGS) -o $@ rapi_wrap.o -L$(dir $(RAPI_BWA)) -lrapi_bwa -L$(BWA_PATH) -lbwa

bwa:
	@echo "BWA_PATH is $(BWA_PATH)"
	$(if $(BWA_PATH),, $(error "You need to set the BWA_PATH variable on the cmd line to point to the compiled BWA source code (e.g., make BWA_PATH=/tmp/bwa)"))


clean:
	rm -f $(SHARED) rapi_wrap.[co] rapi.py rapi.pyc

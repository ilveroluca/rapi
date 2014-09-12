
CC = gcc
AR = ar


PYTHON_INCLUDES := $(shell python-config --includes)
PYTHON_PREFIX := $(shell python-config --prefix)
PYTHON_LIBS := $(shell python-config --libs)

INCLUDES := -I../include $(PYTHON_INCLUDES)
CFLAGS := -g -Wall -O2 -std=c99 -fPIC
RAPI_BWA := ../rapi_bwa/librapi_bwa.a
SHARED := _rapi.so
UNAME := $(shell uname)


.SUFFIXES:.c .o
.PHONY: clean

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDES) $(DFLAGS) $< -o $@

all: $(SHARED)

rapi_wrap.c: ../include/rapi.h rapi.i
	swig -builtin -python $(INCLUDES) -o rapi_wrap.c rapi.i

rapi_wrap.o: rapi_wrap.c
	$(CC) $(CFLAGS) $(INCLUDES) -c rapi_wrap.c

$(SHARED): rapi_wrap.o
# to build the shared library we link against the static rapi_lib
ifeq ($(UNAME), Darwin)
	$(CC) -shared $(CFLAGS) -o $@ rapi_wrap.o -L$(dir $(RAPI_BWA)) -lrapi_bwa -L$(PYTHON_PREFIX)/lib $(PYTHON_LIBS) -lz
	install_name_tool -change libpython2.7.dylib $(PYTHON_PREFIX)/lib/libpython2.7.dylib $(SHARED)
else
	$(CC) -shared $(CFLAGS) -o $@ rapi_wrap.o -L$(dir $(RAPI_BWA)) -lrapi_bwa -lz
endif

clean:
	rm -f $(SHARED) rapi_wrap.[co] rapi.py rapi.pyc

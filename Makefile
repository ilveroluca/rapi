
# Top-level Makefile for rapi project.


# Exported variables are passed down to recursive make calls
# XXX: if you change the value of BWA_PATH change it in the "clean"
# rule as well.
export BWA_PATH := $(PWD)/bwa-auto-build
#
# You can override the value of BWA_PATH by specifying on the command line, like
#     make BWA_PATH=${PWD}/my_other_bwa_dir

$(info "Using BWA_PATH = $(BWA_PATH)")

all: rapi_bwa pyrapi example

bwa_lib: $(BWA_PATH)/libbwa.a

$(BWA_PATH)/libbwa.a:
	rapi_bwa/setup_bwa.sh $(BWA_PATH)

rapi_bwa: bwa_lib
	$(MAKE) -C rapi_bwa/
   
pyrapi: bwa_lib rapi_bwa
	(cd pyrapi && python setup.py build)

example: pyrapi
	$(MAKE) -C example

clean:
	$(MAKE) -C rapi_bwa/ clean
	(cd pyrapi && python setup.py clean --all)


distclean: clean
	# Remove automatically built BWA, if it exists
	rm -rf "$(PWD)/bwa-auto-build"

tests: pyrapi
	python tests/pyrapi/test_pyrapi.py


.PHONY: clean distclean tests pyrapi rapi_bwa example


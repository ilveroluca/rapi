

all: pyrapi


rapi_bwa:
	make -C rapi_bwa/
   
pyrapi: rapi_bwa
	make -C pyrapi/

clean:
	make -C rapi_bwa/ clean
	make -C pyrapi/ clean


tests: pyrapi
	python tests/pyrapi/test_pyrapi.py


.PHONY: clean tests pyrapi rapi_bwa


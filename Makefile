
# Top-level Makefile for rapi project.

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



# Exported variables are passed down to recursive make calls
# XXX: if you change the value of BWA_PATH change it in the "clean"
# rule as well.
export BWA_PATH := $(PWD)/bwa-auto-build
#
# You can override the value of BWA_PATH by specifying on the command line, like
#     make BWA_PATH=${PWD}/my_other_bwa_dir

$(info "Using BWA_PATH = $(BWA_PATH)")

all: rapi_bwa pyrapi jrapi example

bwa_lib: $(BWA_PATH)/libbwa.a

$(BWA_PATH)/libbwa.a:
	rapi_bwa/setup_bwa.sh $(BWA_PATH)

rapi_bwa: bwa_lib
	$(MAKE) -C rapi_bwa/
   
pyrapi: bwa_lib rapi_bwa
	(cd bindings/pyrapi && python setup.py clean --all && python setup.py build)

jrapi: bwa_lib rapi_bwa
	make -C bindings/jrapi

example: pyrapi
	$(MAKE) -C example

clean:
	$(MAKE) -C rapi_bwa/ clean
	$(MAKE) -C bindings/ clean

distclean: clean
	# Remove automatically built BWA, if it exists
	rm -rf "$(PWD)/bwa-auto-build"

tests: pyrapi
	python tests/pyrapi/test_pyrapi.py

.PHONY: clean distclean tests pyrapi rapi_bwa example


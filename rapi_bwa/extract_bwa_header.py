#!/usr/bin/env python

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

import re
import os
import sys

RequiredPrototypes = {
    'kt_for': 'kthread.c',
    'mem_align1_core': 'bwamem.c',
    'mem_approx_mapq_se': 'bwamem.c',
    'mem_mark_primary_se': 'bwamem.c',
    'mem_matesw': 'bwamem_pair.c',
    'mem_pair': 'bwamem_pair.c'
}

def writeline(txt=''):
    sys.stdout.write(txt + '\n')

def extract_prototype(fn_name, filename):
    if os.path.getsize(filename) > 1000000:
        raise RuntimeError("file %s is larger than 1 MB.  Are you sure it's a source file?" % filename)

    with open(filename) as f:
        text = f.read()
    m = re.search(r'([A-Za-z_]\w* %s\([^;{]+)\s*(;|{|//)' % fn_name, text, re.MULTILINE)
    if m:
        proto = m.group(1)
        # remove any trailing comments
        comment_pos = proto.find('//')
        if comment_pos >= 0:
            proto = proto[0:comment_pos]
        # strip whitespace
        proto = proto.strip()
        return proto
    else:
        raise RuntimeError("Couldn't extract prototype for function %s from file %s" % \
                (fn_name, filename))

def extract_bwa_version(bwa_src_path):
    main_file = os.path.join(bwa_src_path, 'main.c')
    with open(main_file) as f:
        for line in f:
            m = re.match(r'#define PACKAGE_VERSION (.+)', line)
            if m:
                bwa_ver = m.group(1)
                return bwa_ver
    raise RuntimeError("Couldn't find PACKAGE_VERSION define in main.c")

def main(args=None):
    if args is None:
        args = sys.argv[1:]

    if args and len(args) != 1:
        sys.exit("Usage: %s [BWA_PATH]" % os.path.basename(sys.argv[0]))

    bwa_path = args[0] if args else '.'

    writeline('#ifndef RAPI_BWA_HEADER')
    writeline('#define RAPI_BWA_HEADER')
    writeline()

    bwa_ver = extract_bwa_version(bwa_path)
    writeline("#define WRAPPED_BWA_VERSION %s" % bwa_ver)
    writeline()

    for fn_name, file_name in RequiredPrototypes.iteritems():
        proto = extract_prototype(fn_name, os.path.join(bwa_path, file_name))
        writeline("extern %s;" % proto)

    writeline()
    writeline('#endif')

if __name__ == '__main__':
    main()

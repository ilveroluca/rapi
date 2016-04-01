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

import sys

def report(*msg):
    print >> sys.stderr, ' '.join(map(str, msg))

class SamFile(object):
    def __init__(self, fp):
        self.fp = fp
        self._lineno = 0
        self._buf = None
        self.advance()

    def advance(self):
        self._buf = self.fp.readline()
        if self._buf:
            self._buf = self._buf.rstrip('\n')
            self._lineno += 1
        else:
            self._buf = None
        return self._buf is not None

    @property
    def current_line(self):
        return self._buf

    @property
    def done(self):
        return self._buf is None

    def read_header(self):
        header_lines = []
        while not self.done and self.current_line.startswith("@"):
            items = self.current_line[1:].split('\t')
            header_lines.append(dict(
                type=items[0],
                tags=items[1:]
                ))
            self.advance()
        return header_lines

def compare_headers(header_a, header_b):
    a = [ h for h in header_a if h['type'] in ('HD', 'SQ') ]
    b = [ h for h in header_b if h['type'] in ('HD', 'SQ') ]
    if a != b:
        if len(a) != len(b):
            report("lengths of headers are different (%d vs. %d)" % (len(a), len(b)))
        else:
            for idx in xrange(len(a)):
                if a[idx] != b[idx]:
                    report("Headers lines are different\n\t%s\nvs\n\t%s\n" % (a[idx], b[idx]))
        return False
    else:
        return True

def compare_alignments(aln1, aln2):
    # expect everything to be the same except tag order
    def find_nth(s, n, char):
        if n <= 0 or not char:
            raise ValueError()
        pos = 0
        while True:
            pos = s.find(char, pos)
            n = n - 1
            if n == 0:
                return pos
            if pos < 0: # not found
                return -1
            pos += 1

    def split_parts(aln):
        start_tags = find_nth(aln, 11, "\t")
        if start_tags < 0: # no tags
            start_tags = len(aln)
        tags = set( t for t in aln[start_tags+1:].split('\t') )
        return dict(aln=aln[0:start_tags], tags=tags)

    alignments = map(split_parts, (aln1, aln2))

    result = False
    if alignments[0]['aln'] != alignments[1]['aln']:
        report("alignments DIFFER: \n\t'%s'\nand\n\t'%s'" % (alignments[0]['aln'], alignments[1]['aln']))
    elif alignments[0]['tags'] != alignments[1]['tags']:
        aln_ids = [ a['aln'].split('\t', 1)[0] for a in alignments ]
        report("alignment tags DIFFER for alignments %s and %s -> { %s } != { %s }" %
                (aln_ids[0], aln_ids[1], ','.join(alignments[0]['tags']), ','.join(alignments[1]['tags'])))
    else:
        result = True
    return result

def compare_sam_files(file_a, file_b):
    with open(file_a) as a, open(file_b) as b:
        sams = map(SamFile, (a, b))
        headers = [ v.read_header() for v in sams ]
        if not compare_headers(*headers):
           report("headers are different")
           return False

        all_ok = True
        while not any(s.done for s in sams):
            this_ok = compare_alignments(sams[0].current_line, sams[1].current_line)
            all_ok = all_ok and this_ok
            for s in sams:
                s.advance()

        if not all(s.done for s in sams):
            report("Finished reading one SAM before the other!")
            report("file_a is %s done; file_b is %s done" %
                    ( '' if sams[0].done else 'NOT', '' if sams[1].done else 'NOT'))
            all_ok = False
    return all_ok

def main(args=None):
    if args is None:
        args = sys.argv[1:]

    file_a, file_b = args
    ok = compare_sam_files(file_a, file_b)
    sys.exit(0 if ok else 1)

if __name__ == '__main__':
    main()

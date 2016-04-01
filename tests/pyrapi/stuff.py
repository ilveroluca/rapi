

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

import os

MiniRef = \
        os.path.abspath(
            os.path.join(
                os.path.dirname(__file__), '../mini_ref/mini_ref.fasta'))

MiniRefSequencesTxt = \
        os.path.abspath(
            os.path.join(
                os.path.dirname(__file__), '../mini_ref/mini_ref_seqs.txt'))

# we cache the list produced by get_mini_ref_seqs
_mini_ref_seqs = None

def read_seqs(filename):
    """
    Fetch a tab-separated file of sequences::

        ID	SEQ1	Q1	SEQ2	Q2
    """
    with open(filename) as f:
        # use a tuple since it's immutable
        seqs = tuple([ line.rstrip('\n').split('\t') for line in f ])
    return seqs

def get_mini_ref_seqs():
    global _mini_ref_seqs
    if _mini_ref_seqs is None:
        _mini_ref_seqs = read_seqs(MiniRefSequencesTxt)
    return _mini_ref_seqs


_complement = {
    'A':'T',
    'G':'C',
    'C':'G',
    'T':'A',
    'N':'N',
    'a':'t',
    'g':'c',
    'c':'g',
    't':'a',
    'n':'n'
}

def rev_complement(seq):
    return ''.join( _complement[base] for base in seq[::-1] )

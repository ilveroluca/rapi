

import os
import sys

src_root_dir = \
        os.path.abspath(
            os.path.join(
                os.path.dirname(__file__), '../..'))

def append_src_to_pypath():
    sys.path.append(src_root_dir)


MiniRef = \
        os.path.abspath(
            os.path.join(
                os.path.dirname(__file__), '../mini_ref/mini_ref.fasta'))

MiniRefSequencesTxt = \
        os.path.abspath(
            os.path.join(
                os.path.dirname(__file__), '../mini_ref/mini_ref_seqs.txt'))


SequencesTxt = \
        os.path.abspath(
            os.path.join(
                os.path.dirname(__file__), '../seqs.txt'))

# we cache the list produced by get_sequences
_mini_ref_seqs = None
_seqs = None

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

def get_sequences():
    global _seqs
    if _seqs is None:
        _seqs = read_seqs(SequencesTxt)
    return _seqs



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

SequencesTxt = \
        os.path.abspath(
            os.path.join(
                os.path.dirname(__file__), '../seqs.txt'))

# we cache the list produced by get_sequences
_seqs = None

def get_sequences():
    """
    Fetch a tab-separated file of sequences::

        ID	SEQ1	Q1	SEQ2	Q2
    """
    global _seqs
    if _seqs is None:
        with open(SequencesTxt) as f:
            # use a tuple since it's immutable
            _seqs = tuple([ line.rstrip('\n').split('\t') for line in f ])
    return _seqs

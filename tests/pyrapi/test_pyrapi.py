#!/usr/bin/env python

import os
import sys
import unittest

sys.path.append(
        os.path.abspath(
            os.path.join(
                os.path.dirname(__file__), '../..')))

import pyrapi.rapi as rapi


MiniRef = \
        os.path.abspath(
            os.path.join(
                os.path.dirname(__file__), '../mini_ref/mini_ref.fasta'))

SequencesTxt = \
        os.path.abspath(
            os.path.join(
                os.path.dirname(__file__), '../seqs.txt'))

def read_sequences():
    """
    Read a tab-separated file of sequences::

        ID	SEQ1	Q1	SEQ2	Q2
    """
    with open(SequencesTxt) as f:
        return [ line.rstrip('\n').split('\t') for line in f ]

class TestPyrapi(unittest.TestCase):
    def setUp(self):
        self.opts = rapi.opts()
        rapi.init(self.opts)

    def tearDown(self):
        rapi.shutdown()

    def test_init(self):
        r = rapi.init(self.opts)
        self.assertIsNone(r)

    def test_shutdown(self):
        r = rapi.shutdown()
        self.assertIsNone(r)

    def test_bad_ref_path(self):
        self.assertRaises(RuntimeError, rapi.ref, 'bad/path')

    def test_ref_load_unload(self):
        r = rapi.ref(MiniRef)
        self.assertEqual(MiniRef, r.path)
        self.assertEqual(1, len(r))
        r.unload()

    def test_ref_contigs(self):
        r = rapi.ref(MiniRef)
        try:
            contigs = [ c for c in r ]
            self.assertEquals(1, len(contigs))
            self.assertEquals(len(r), len(contigs))

            c0 = contigs[0]
            self.assertEquals('chr1', c0.name)
            self.assertEquals(60000, c0.len)
            for v in (c0.assembly_identifier, c0.species, c0.uri, c0.md5):
                self.assertIsNone(v)
        finally:
            r.unload()

    def test_batch_wrap_create_(self):
        w = rapi.batch_wrap(2)
        self.assertEquals(2, w.n_reads_per_frag())
        self.assertEquals(0, len(w))

        w = rapi.batch_wrap(1)
        self.assertEquals(1, w.n_reads_per_frag())

    def test_batch_wrap_reserve(self):
        w = rapi.batch_wrap(2)
        w.reserve(10)
        self.assertGreaterEqual(10, w.capacity())
        w.reserve(100)
        self.assertGreaterEqual(100, w.capacity())

    def test_batch_wrap_append_one(self):
        w = rapi.batch_wrap(2)
        seq_pair = read_sequences()[0]
        # insert one read
        w.append(seq_pair[0], seq_pair[1], seq_pair[2], rapi.QENC_SANGER)
        self.assertEquals(1, len(w))
        read1 = w.get_read(0, 0)
        self.assertEquals(seq_pair[0], read1.id)
        self.assertEquals(seq_pair[1], read1.seq)
        self.assertEquals(seq_pair[2], read1.qual)
        self.assertEquals(len(seq_pair[1]), len(read1))

    def test_batch_wrap_append_illumina(self):
        w = rapi.batch_wrap(2)
        seq_pair = read_sequences()[0]
        # convert the base qualities from sanger to illumina encoding
        # (subtract sanger offset and add illumina offset)
        illumina_qualities = ''.join([ chr(ord(c) - 33 + 64) for c in seq_pair[2] ])
        w.append(seq_pair[0], seq_pair[1], illumina_qualities, rapi.QENC_ILLUMINA)
        read1 = w.get_read(0, 0)
        # Internally base qualities should be converted to sanger format
        self.assertEquals(seq_pair[2], read1.qual)

    def test_batch_wrap_append_baseq_out_of_range(self):
        w = rapi.batch_wrap(2)
        seq_pair = read_sequences()[0]

        new_q = chr(32)*len(seq_pair[1]) # sanger encoding goes down to 33
        self.assertRaises(ValueError, w.append, seq_pair[0], seq_pair[1], new_q, rapi.QENC_SANGER)

        new_q = chr(127)*len(seq_pair[1]) # and up to 126
        self.assertRaises(ValueError, w.append, seq_pair[0], seq_pair[1], new_q, rapi.QENC_SANGER)

        new_q = chr(63)*len(seq_pair[1]) # illumina encoding goes down to 64
        self.assertRaises(ValueError, w.append, seq_pair[0], seq_pair[1], new_q, rapi.QENC_ILLUMINA)

    def test_n_reads_per_frag(self):
        w = rapi.batch_wrap(2)
        self.assertEquals(2, w.n_reads_per_frag())
        w = rapi.batch_wrap(1)
        self.assertEquals(1, w.n_reads_per_frag())


def suite():
    return unittest.TestLoader().loadTestsFromTestCase(TestPyrapi)

def main():
    result = unittest.TextTestRunner(verbosity=2).run(suite())
    if not result.wasSuccessful():
        sys.exit(1)

if __name__ == '__main__':
    main()

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
        r = rapi.init(rapi.opts())
        self.assertIsNone(r)

    def test_shutdown(self):
        rapi.init(rapi.opts())
        r = rapi.shutdown()
        self.assertIsNone(r)

    def test_aligner_name_vers(self):
        name = rapi.aligner_name()
        ver = rapi.aligner_version()
        self.assertGreater(len(name), 0)
        self.assertGreater(len(ver), 0)


class TestPyrapiRef(unittest.TestCase):
    def setUp(self):
        self.opts = rapi.opts()
        rapi.init(self.opts)
        self.ref = rapi.ref(MiniRef)

    def tearDown(self):
        self.ref.unload()
        rapi.shutdown()

    def test_bad_ref_path(self):
        self.assertRaises(RuntimeError, rapi.ref, 'bad/path')

    def test_ref_load_unload(self):
        self.assertEqual(MiniRef, self.ref.path)
        self.assertEqual(1, len(self.ref))

    def test_ref_len(self):
        self.assertEqual(1, len(self.ref))

    def test_ref_iteration(self):
        contigs = [ c for c in self.ref ]
        self.assertEquals(1, len(contigs))
        self.assertEquals(len(self.ref), len(contigs))

        c0 = contigs[0]
        self.assertEquals('chr1', c0.name)
        self.assertEquals(60000, c0.len)
        for v in (c0.assembly_identifier, c0.species, c0.uri, c0.md5):
            self.assertIsNone(v)

    def test_get_contig(self):
        c0 = self.ref.get_contig(0)
        self.assertEquals('chr1', c0.name)
        self.assertEquals(60000, c0.len)

    def test_get_contig_out_of_bounds(self):
        self.assertRaises(IndexError, self.ref.get_contig, 1)
        self.assertRaises(IndexError, self.ref.get_contig, -1)

    def test_get_contig_by_get_item(self):
        c0 = self.ref[0]
        self.assertEquals('chr1', c0.name)
        self.assertEquals(60000, c0.len)

    def test_get_contig_by_get_item_out_of_bounds(self):
        self.assertRaises(IndexError, self.ref.__getitem__, 1)
        self.assertRaises(IndexError, self.ref.__getitem__, -1)


class TestPyrapiReadBatch(unittest.TestCase):
    def setUp(self):
        self.opts = rapi.opts()
        rapi.init(self.opts)
        self.w = rapi.read_batch(2)

    def tearDown(self):
        rapi.shutdown()

    def test_create_bad_arg(self):
        self.assertRaises(ValueError, rapi.read_batch, -1)

    def test_create(self):
        self.assertEquals(2, self.w.n_reads_per_frag())
        self.assertEquals(0, len(self.w))

        w = rapi.read_batch(1)
        self.assertEquals(1, w.n_reads_per_frag())

    def test_reserve(self):
        self.w.reserve(10)
        self.assertGreaterEqual(10, self.w.capacity())
        self.w.reserve(100)
        self.assertGreaterEqual(100, self.w.capacity())
        self.assertRaises(ValueError, self.w.reserve, -1)

    def test_append_one(self):
        seq_pair = read_sequences()[0]
        # insert one read
        self.w.append(seq_pair[0], seq_pair[1], seq_pair[2], rapi.QENC_SANGER)
        self.assertEquals(1, len(self.w))
        read1 = self.w.get_read(0, 0)
        self.assertEquals(seq_pair[0], read1.id)
        self.assertEquals(seq_pair[1], read1.seq)
        self.assertEquals(seq_pair[2], read1.qual)
        self.assertEquals(len(seq_pair[1]), len(read1))

    def test_append_illumina(self):
        seq_pair = read_sequences()[0]
        # convert the base qualities from sanger to illumina encoding
        # (subtract sanger offset and add illumina offset)
        illumina_qualities = ''.join([ chr(ord(c) - 33 + 64) for c in seq_pair[2] ])
        self.w.append(seq_pair[0], seq_pair[1], illumina_qualities, rapi.QENC_ILLUMINA)
        read1 = self.w.get_read(0, 0)
        # Internally base qualities should be converted to sanger format
        self.assertEquals(seq_pair[2], read1.qual)

    def test_append_baseq_out_of_range(self):
        seq_pair = read_sequences()[0]

        new_q = chr(32)*len(seq_pair[1]) # sanger encoding goes down to 33
        self.assertRaises(ValueError, self.w.append, seq_pair[0], seq_pair[1], new_q, rapi.QENC_SANGER)

        new_q = chr(127)*len(seq_pair[1]) # and up to 126
        self.assertRaises(ValueError, self.w.append, seq_pair[0], seq_pair[1], new_q, rapi.QENC_SANGER)

        new_q = chr(63)*len(seq_pair[1]) # illumina encoding goes down to 64
        self.assertRaises(ValueError, self.w.append, seq_pair[0], seq_pair[1], new_q, rapi.QENC_ILLUMINA)

    def test_n_reads_per_frag(self):
        self.assertEquals(2, self.w.n_reads_per_frag())
        w = rapi.read_batch(1)
        self.assertEquals(1, w.n_reads_per_frag())

    def test_get_read_out_of_bounds(self):
        self.assertEqual(0, len(self.w))
        self.assertRaises(IndexError, self.w.get_read, 0, 0)
        self.assertRaises(IndexError, self.w.get_read, -1, 0)
        self.assertRaises(IndexError, self.w.get_read, -1, -1)

        # now load a sequence and ensure we get exceptions if we access beyond the limits
        seq_pair = read_sequences()[0]
        self.w.append(seq_pair[0], seq_pair[1], seq_pair[2], rapi.QENC_SANGER)
        self.assertIsNotNone(self.w.get_read(0, 0))
        self.assertRaises(IndexError, self.w.get_read, 0, 1)
        self.assertRaises(IndexError, self.w.get_read, 1, 0)
        self.assertRaises(IndexError, self.w.get_read, 1, 1)




def suite():
    s = unittest.TestLoader().loadTestsFromTestCase(TestPyrapi)
    s.addTests(unittest.TestLoader().loadTestsFromTestCase(TestPyrapiRef))
    s.addTests(unittest.TestLoader().loadTestsFromTestCase(TestPyrapiReadBatch))
    return s

def main():
    result = unittest.TextTestRunner(verbosity=2).run(suite())
    if not result.wasSuccessful():
        sys.exit(1)

if __name__ == '__main__':
    main()

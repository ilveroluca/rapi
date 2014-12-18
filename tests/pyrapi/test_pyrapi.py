#!/usr/bin/env python

import re
import sys
import unittest

import stuff
stuff.append_src_to_pypath()

import pyrapi.rapi as rapi

class TestPyrapi(unittest.TestCase):
    def setUp(self):
        self.opts = rapi.opts()
        rapi.init(self.opts)

    def tearDown(self):
        rapi.shutdown()

    def test_init(self):
        r = rapi.init(rapi.opts())
        self.assertIsNone(r)

    def test_init_opts_none(self):
        rapi.init(None)

    def test_shutdown(self):
        rapi.init(rapi.opts())
        r = rapi.shutdown()
        self.assertIsNone(r)

    def test_aligner_name_vers(self):
        name = rapi.aligner_name()
        ver = rapi.aligner_version()
        plugin_version = rapi.plugin_version()
        self.assertGreater(len(name), 0)
        self.assertGreater(len(ver), 0)
        self.assertGreater(len(plugin_version), 0)


class TestPyrapiRef(unittest.TestCase):
    def setUp(self):
        self.opts = rapi.opts()
        rapi.init(self.opts)
        self.ref = rapi.ref(stuff.MiniRef)

    def tearDown(self):
        self.ref.unload()
        rapi.shutdown()

    def test_bad_ref_path(self):
        self.assertRaises(RuntimeError, rapi.ref, 'bad/path')
        self.assertRaises(TypeError, rapi.ref, None)

    def test_ref_load_unload(self):
        self.assertEqual(stuff.MiniRef, self.ref.path)
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
        self.assertRaises(TypeError, self.ref.get_contig, None)


    def test_get_contig_by_get_item(self):
        c0 = self.ref[0]
        self.assertEquals('chr1', c0.name)
        self.assertEquals(60000, c0.len)

    def test_get_contig_by_get_item_out_of_bounds(self):
        self.assertRaises(IndexError, self.ref.__getitem__, 1)
        self.assertRaises(IndexError, self.ref.__getitem__, -1)
        self.assertRaises(TypeError, self.ref.__getitem__, None)


class TestPyrapiReadBatch(unittest.TestCase):
    def setUp(self):
        self.opts = rapi.opts()
        rapi.init(self.opts)
        self.w = rapi.read_batch(2)

    def tearDown(self):
        rapi.shutdown()

    def test_hidden_members(self):
        self.assertFalse(hasattr(self.w, 'batch'))
        self.assertFalse(hasattr(self.w, 'fragment_num'))

    def test_create_bad_arg(self):
        self.assertRaises(ValueError, rapi.read_batch, -1)
        self.assertRaises(TypeError, rapi.read_batch, None)

    def test_create(self):
        self.assertEquals(2, self.w.n_reads_per_frag)
        self.assertEquals(0, len(self.w))
        self.assertEquals(0, self.w.n_fragments)

        w = rapi.read_batch(1)
        self.assertEquals(1, w.n_reads_per_frag)

    def test_reserve(self):
        self.w.reserve(10)
        self.assertGreaterEqual(10, self.w.capacity)
        self.w.reserve(100)
        self.assertGreaterEqual(100, self.w.capacity)
        self.assertEquals(0, self.w.n_fragments)
        self.assertEquals(0, len(self.w))
        self.assertRaises(ValueError, self.w.reserve, -1)
        self.assertRaises(TypeError, self.w.reserve, None)

    def test_append_one(self):
        seq_pair = stuff.get_mini_ref_seqs()[0]
        # insert one read
        self.w.append(seq_pair[0], seq_pair[1], seq_pair[2], rapi.QENC_SANGER)
        self.assertEquals(1, len(self.w))
        read1 = self.w.get_read(0, 0)
        self.assertEquals(seq_pair[0], read1.id)
        self.assertEquals(seq_pair[1], read1.seq)
        self.assertEquals(seq_pair[2], read1.qual)
        self.assertEquals(len(seq_pair[1]), len(read1))
        self.assertEquals(0, self.w.n_fragments)

    def test_append_without_qual(self):
        seq_pair = stuff.get_mini_ref_seqs()[0]
        self.w.append(seq_pair[0], seq_pair[1], None, rapi.QENC_SANGER)
        self.assertEquals(1, len(self.w))
        read1 = self.w.get_read(0, 0)
        self.assertEquals(seq_pair[1], read1.seq)
        self.assertIsNone(read1.qual)

    def test_append_illumina(self):
        seq_pair = stuff.get_mini_ref_seqs()[0]
        # convert the base qualities from sanger to illumina encoding
        # (subtract sanger offset and add illumina offset)
        illumina_qualities = ''.join([ chr(ord(c) - 33 + 64) for c in seq_pair[2] ])
        self.w.append(seq_pair[0], seq_pair[1], illumina_qualities, rapi.QENC_ILLUMINA)
        read1 = self.w.get_read(0, 0)
        # Internally base qualities should be converted to sanger format
        self.assertEquals(seq_pair[2], read1.qual)

    def test_append_baseq_out_of_range(self):
        seq_pair = stuff.get_mini_ref_seqs()[0]

        new_q = chr(32)*len(seq_pair[1]) # sanger encoding goes down to 33
        self.assertRaises(ValueError, self.w.append, seq_pair[0], seq_pair[1], new_q, rapi.QENC_SANGER)

        new_q = chr(127)*len(seq_pair[1]) # and up to 126
        self.assertRaises(ValueError, self.w.append, seq_pair[0], seq_pair[1], new_q, rapi.QENC_SANGER)

        new_q = chr(63)*len(seq_pair[1]) # illumina encoding goes down to 64
        self.assertRaises(ValueError, self.w.append, seq_pair[0], seq_pair[1], new_q, rapi.QENC_ILLUMINA)

    def test_append_none(self):
        self.assertRaises(ValueError, self.w.append, None, None, None, rapi.QENC_SANGER)
        self.assertRaises(ValueError, self.w.append, "some id", None, None, rapi.QENC_SANGER)
        self.assertRaises(TypeError, self.w.append, "some id", "AGCT", None, None)

    def test_clear(self):
        seq_pair = stuff.get_mini_ref_seqs()[0]
        self.w.append(seq_pair[0], seq_pair[1], seq_pair[2], rapi.QENC_SANGER)
        self.assertEquals(1, len(self.w))
        self.w.clear()
        # after clear, the length should be back to zero, but the capacity should be unaffected.
        self.assertEquals(0, len(self.w))
        self.assertGreaterEqual(self.w.capacity, 1)

    def test_set_read_out_of_bounds(self):
        self.assertRaises(ValueError, self.w.set_read, 0, 0, "some id", "AGCT", None, rapi.QENC_SANGER)
        self.w.reserve(2) # 2 reads, 1 fragment
        self.w.set_read(0, 1, "some id", "AGCT", None, rapi.QENC_SANGER)
        self.assertEqual(2, len(self.w))
        self.assertRaises(ValueError, self.w.set_read, 1, 0, "some id", "AGCT", None, rapi.QENC_SANGER)
        self.assertRaises(ValueError, self.w.set_read, 1, 1, "some id", "AGCT", None, rapi.QENC_SANGER)

    def test_set_read_none(self):
        self.w.reserve(1)
        self.assertRaises(ValueError, self.w.set_read, 0, 0, None, None, None, rapi.QENC_SANGER)
        self.assertRaises(ValueError, self.w.set_read, 0, 0, "some id", None, None, rapi.QENC_SANGER)
        self.assertRaises(TypeError, self.w.set_read, 0, 0, "some id", "AGCT", None, None)
        self.assertRaises(TypeError, self.w.set_read, None, 0, "some id", "AGCT", None, rapi.QENC_SANGER)
        self.assertRaises(TypeError, self.w.set_read, 0, None, "some id", "AGCT", None, rapi.QENC_SANGER)

    def test_n_reads_per_frag(self):
        self.assertEquals(2, self.w.n_reads_per_frag)
        w = rapi.read_batch(1)
        self.assertEquals(1, w.n_reads_per_frag)

    def test_n_fragments(self):
        seq_pair = stuff.get_mini_ref_seqs()[0]
        self.assertEquals(0, self.w.n_fragments)
        # insert one read
        self.w.append(seq_pair[0], seq_pair[1], seq_pair[2], rapi.QENC_SANGER)
        self.assertEquals(0, self.w.n_fragments)
        self.w.append(seq_pair[0], seq_pair[3], seq_pair[4], rapi.QENC_SANGER)
        self.assertEquals(1, self.w.n_fragments)

    def test_get_read_out_of_bounds(self):
        self.assertEqual(0, len(self.w))
        self.assertRaises(IndexError, self.w.get_read, 0, 0)
        self.assertRaises(IndexError, self.w.get_read, -1, 0)
        self.assertRaises(IndexError, self.w.get_read, -1, -1)

        # now load a sequence and ensure we get exceptions if we access beyond the limits
        seq_pair = stuff.get_mini_ref_seqs()[0]
        self.w.append(seq_pair[0], seq_pair[1], seq_pair[2], rapi.QENC_SANGER)
        self.assertIsNotNone(self.w.get_read(0, 0))
        self.assertRaises(IndexError, self.w.get_read, 0, 1)
        self.assertRaises(IndexError, self.w.get_read, 1, 0)
        self.assertRaises(IndexError, self.w.get_read, 1, 1)
        self.assertRaises(TypeError, self.w.get_read, None, 0)
        self.assertRaises(TypeError, self.w.get_read, 0, None)

    def test_iteration_values(self):
        seqs = stuff.get_mini_ref_seqs()[0:2]
        for pair in seqs:
            self.w.append(pair[0], pair[1], pair[2], rapi.QENC_SANGER)
            self.w.append(pair[0], pair[3], pair[4], rapi.QENC_SANGER)

        for idx, fragment in enumerate(self.w):
            self.assertEquals(2, len(fragment))
            self.assertEquals(seqs[idx][0], fragment[0].id)
            self.assertEquals(seqs[idx][1], fragment[0].seq)
            self.assertEquals(seqs[idx][0], fragment[1].id)
            self.assertEquals(seqs[idx][3], fragment[1].seq)

    def test_iteration_completeness(self):
        # empty batch
        self.assertEquals(0, sum(1 for frag in self.w))
        # now put some sequences into it
        seqs = stuff.get_mini_ref_seqs()
        for pair in seqs:
            self.w.append(pair[0], pair[1], pair[2], rapi.QENC_SANGER)
            self.w.append(pair[0], pair[3], pair[4], rapi.QENC_SANGER)
        # test whether we iterate over all the fragments
        self.assertEquals(len(seqs), sum(1 for frag in self.w))

    def test_iteration_completeness_single_end(self):
        seqs = stuff.get_mini_ref_seqs()
        batch = rapi.read_batch(1)
        for pair in seqs:
            batch.append(pair[0], pair[1], pair[2], rapi.QENC_SANGER)
            batch.append(pair[0], pair[3], pair[4], rapi.QENC_SANGER)
        # test whether we iterate over all the fragments
        self.assertEquals(2 * len(seqs), sum(1 for frag in batch))
        it = iter(batch)
        fragment = next(it)
        self.assertEquals(1, len(fragment))


class TestPyrapiAlignment(unittest.TestCase):
    def setUp(self):
        self.opts = rapi.opts()
        rapi.init(self.opts)
        self.ref = rapi.ref(stuff.MiniRef)
        aligner = rapi.aligner(self.opts)
        self.batch = rapi.read_batch(2)
        reads = stuff.get_mini_ref_seqs()
        for row in reads:
            self.batch.append(row[0], row[1], row[2], rapi.QENC_SANGER)
            self.batch.append(row[0], row[3], row[4], rapi.QENC_SANGER)
        aligner.align_reads(self.ref, self.batch)

    def tearDown(self):
        self.ref.unload()
        rapi.shutdown()

    def test_alignment_struct(self):
        rapi_read = self.batch.get_read(0, 0)
        self.assertEqual("read_00", rapi_read.id)

        self.assertTrue(rapi_read.n_alignments > 0)

        aln = rapi_read.get_aln(0)
        self.assertEqual("chr1", aln.contig.name)
        self.assertEqual(32461, aln.pos)
        self.assertEqual("60M", aln.get_cigar_string())
        self.assertEqual( ( ('M', 60), ), aln.get_cigar_ops())
        self.assertEqual(60, aln.mapq)
        self.assertEqual(60, aln.score)
        # Flag 65 is 'p1'
        self.assertTrue(aln.paired)
        self.assertFalse(aln.prop_paired)
        self.assertTrue(aln.mapped)
        self.assertFalse(aln.reverse_strand)
        self.assertFalse(aln.secondary_aln)
        self.assertEqual(0, aln.n_mismatches)
        self.assertEqual(0, aln.n_gap_opens)
        self.assertEqual(0, aln.n_gap_extensions)

        expected_tags = dict(
                MD='60',
                XS=0)
        # the SAM has more tags (NM and AS), but RAPI makes that information
        # available through other members of the alignment structure.
        self.assertEqual(expected_tags, aln.get_tags())

        # check out reads that don't align perfectly
        aln = self.batch.get_read(1, 0).get_aln(0)
        self.assertEqual(60, aln.mapq)
        self.assertEqual(51, aln.score)
        self.assertEqual('11M3D49M', aln.get_cigar_string())
        self.assertEqual( (('M', 11), ('D', 3), ('M', 49)), aln.get_cigar_ops())
        self.assertEqual(3, aln.n_mismatches)
        md_tag = aln.get_tags()['MD']
        self.assertEqual('11^CCC49', md_tag)
        self.assertEqual(0, aln.n_gap_opens)
        self.assertEqual(0, aln.n_gap_extensions)

        aln = self.batch.get_read(2, 0).get_aln(0)
        self.assertEqual(60, aln.mapq)
        self.assertEqual(48, aln.score)
        self.assertEqual('13M3I44M', aln.get_cigar_string())
        self.assertEqual( (('M', 13), ('I', 3), ('M', 44)), aln.get_cigar_ops())
        self.assertEqual(3, aln.n_mismatches)
        md_tag = aln.get_tags()['MD']
        self.assertEqual('57', md_tag)

        aln = self.batch.get_read(3, 0).get_aln(0)
        self.assertEqual(60, aln.mapq)
        self.assertEqual(50, aln.score)
        self.assertEqual('60M', aln.get_cigar_string())
        self.assertEqual(2, aln.n_mismatches)
        md_tag = aln.get_tags()['MD']
        self.assertEqual('15T16C27', md_tag)



    def test_get_aln_out_of_bounds(self):
        rapi_read = self.batch.get_read(0, 0)
        self.assertRaises(IndexError, rapi_read.get_aln, -1)
        self.assertRaises(IndexError, rapi_read.get_aln, rapi_read.n_alignments)

    def test_alignment_iter(self):
        rapi_read = self.batch.get_read(0, 0)
        alignments = [ a for a in rapi_read.iter_aln() ]
        self.assertEqual(rapi_read.n_alignments, len(alignments))
        a = alignments[0]
        self.assertEqual("chr1", a.contig.name)
        self.assertEqual(32461, a.pos)

        ialn = rapi_read.iter_aln()
        b = next(ialn)
        self.assertEqual("chr1", b.contig.name)
        self.assertEqual(32461, b.pos)
        # ensure the iterator is read-only
        def assign_to_n_left(x):
            ialn.n_left = x
        self.assertRaises(AttributeError, assign_to_n_left, 33)

    def _compare_sam_records(self, a_sam, b_sam):
        print >> sys.stderr, "comparing \n\n\t", a_sam, "\n\t", b_sam
        # The two outputs may not be identical because the order or the tags isn't defined.

        a_tag_start = re.search(r'\t[A-Z][A-Z]:[A-z]:.*', a_sam).start()
        b_tag_start = re.search(r'\t[A-Z][A-Z]:[A-z]:.*', b_sam).start()

        self.assertEqual(a_sam[0:a_tag_start], b_sam[0:b_tag_start])

        # get the part of the string that contains the tags.  The '+1' is
        # because the regex we used to find the tags began with a tab, which
        # we want to skip.
        a_tags = set( a_sam[a_tag_start+1:].split('\t') )
        b_tags = set( b_sam[b_tag_start+1:].split('\t') )
        self.assertEquals(a_tags, b_tags)

    def test_sam(self):
        # We ran this command line:
        #     bwa mem -p -T 0 -a mini_ref/mini_ref.fasta mini_ref/mini_ref_seqs.fastq
        # First 8 lines == first two reads from fastq file.
        # BWA version: 0.7.8-r455
        # The BWA options only serve to tell it not to filter any reads from the output:
        #    * -p: first query file consists of interleaved paired-end sequences
        #    * -T 0: minimum score to output [30]
        #    * -a output all alignments for SE or unpaired PE
        #
        # We're using no options with RAPI, which is equivalent to that command line.
        # This is the output from BWA, without the SAM header.
        expected = [
            """read_00	65	chr1	32461	60	60M	=	32581	121	AAAACTGACCCACACAGAAAAACTAATTGTGAGAACCAATATTATACTAAATTCATTTGA	EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	NM:i:0	MD:Z:60	AS:i:60	XS:i:0""",
            """read_00	129	chr1	32581	60	60M	=	32461	-121	CAAAAGTTAACCCATATGGAATGCAATGGAGGAAATCAATGACATATCAGATCTAGAAAC	EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	NM:i:0	MD:Z:60	AS:i:60	XS:i:0""",
            """read_00_rev	113	chr1	32461	60	60M	=	32581	121	AAAACTGACCCACACAGAAAAACTAATTGTGAGAACCAATATTATACTAAATTCATTTGA	EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	NM:i:0	MD:Z:60	AS:i:60	XS:i:0""",
            """read_00_rev	177	chr1	32581	60	60M	=	32461	-121	CAAAAGTTAACCCATATGGAATGCAATGGAGGAAATCAATGACATATCAGATCTAGAAAC	EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	NM:i:0	MD:Z:60	AS:i:60	XS:i:0"""
        ]

        # we produce SAM for the first pair in our set and the last one (which is the
        # first pair reversed and complemented)
        rapi_sam = rapi.format_sam(self.batch, 0).split('\n') + rapi.format_sam(self.batch, 4).split('\n')

        for i in xrange(len(rapi_sam)):
            self._compare_sam_records(expected[i], rapi_sam[i])

    def test_get_insert_size(self):
        aln_read = self.batch.get_read(0, 0).get_aln(0)
        aln_mate = self.batch.get_read(0, 1).get_aln(0)
        self.assertEqual(121, rapi.get_insert_size(aln_read, aln_mate))
        self.assertEqual(-121, rapi.get_insert_size(aln_mate, aln_read))
        self.assertEqual(0, rapi.get_insert_size(aln_read, aln_read))


#    def test_align_se(self):
#        aligner = rapi.aligner(self.opts)
#        batch = rapi.read_batch(1)
#        reads = stuff.get_mini_ref_seqs()
#        batch.append(reads[0][0], reads[0][1], reads[0][2], rapi.QENC_SANGER)
#        batch.append(reads[1][0], reads[1][1], reads[1][2], rapi.QENC_SANGER)
#        aligner.align_reads(self.ref, batch)
#
#        rapi_read = batch.get_read(0, 0)
#        self.assertTrue(rapi_read.n_alignments > 0)
#        rapi_read = batch.get_read(1, 0)
#        self.assertTrue(rapi_read.n_alignments > 0)


def suite():
    s = unittest.TestLoader().loadTestsFromTestCase(TestPyrapi)
    s.addTests(unittest.TestLoader().loadTestsFromTestCase(TestPyrapiRef))
    s.addTests(unittest.TestLoader().loadTestsFromTestCase(TestPyrapiReadBatch))
    s.addTests(unittest.TestLoader().loadTestsFromTestCase(TestPyrapiAlignment))
    return s

def main():
    result = unittest.TextTestRunner(verbosity=2).run(suite())
    return 0 if result.wasSuccessful() else 1

if __name__ == '__main__':
    sys.exit(main())

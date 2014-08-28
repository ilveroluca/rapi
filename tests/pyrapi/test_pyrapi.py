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

print "MiniRef: ", MiniRef


class TestPyrapi(unittest.TestCase):
    def setUp(self):
        self.opts = rapi.opts()
        rapi.init(self.opts)

    def tearDown(self):
        rapi.shutdown()

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

def suite():
    return unittest.TestLoader().loadTestsFromTestCase(TestPyrapi)

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())

#!/usr/bin/env python

import argparse
import logging
import re
import sys

import pyrapi

logging.basicConfig(level=logging.INFO)
_log = logging.getLogger('align')

def _read_fq_record(fp):
    # read 4 lines
    lines = list()
    lineno = 0
    for _ in xrange(4):
        lineno += 1
        line = fp.readline()
        if line:
            lines.append(line.rstrip('\n'))
        else:
            return None
    if lines[2] != '+':
        raise RuntimeError("Format error at line %d:\n%s" % (lineno, '\n'.join(lines)))

    retval = dict()
    retval['id'] = re.sub(r'^@\s*', '', lines[0])
    retval['seq'] = lines[1]
    retval['q'] = lines[3]
    return retval

def read_fastq_1(fp):
    done = False
    try:
        while not done:
            record1 = _read_fq_record(fp)
            record2 = _read_fq_record(fp)
            if record1 and record2:
                yield record1, record2
            else:
                done = True
    finally:
        fp.close()

def read_fastq_2(f1, f2):
    done = False
    try:
        while not done:
            record1 = _read_fq_record(f1)
            record2 = _read_fq_record(f2)
            if record1 and record2:
                yield record1, record2
            else:
                done = True
    finally:
        f1.close()
        f2.close()

def read_prq(fp):
    try:
        for line in fp:
            row = line.rstrip('\n').split('\t')
            record1 = dict( id=row[0], seq=row[1], q=row[2] )
            record2 = dict( id=row[0], seq=row[3], q=row[4] )
            yield record1, record2
    finally:
        fp.close()


def create_input(options):
    if options.format == 'fastq':
        if len(options.input) == 2:
            return read_fastq_2(options.input[0], options.input[1])
        elif len(options.input) == 1:
            return read_fastq_1(options.input[0])
        else:
            raise RuntimeError("Unexpected number of fastq input files %d!" % len(options.input))
    elif options.format == 'prq':
        if len(options.input) == 1:
            return read_prq(options.input[0])
        else:
            raise RuntimeError("Unexpected number of prq input files %d!" % len(options.input))



def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('ref')
    parser.add_argument('input', type=argparse.FileType('r'), nargs='*', default=sys.stdin)
    parser.add_argument('--format', choices=['fastq', 'prq'], default='fastq')
    parser.add_argument('--se', action="store_true", default=False)

    options = parser.parse_args(args)

    if len(options.input) > 2:
        parser.error("Number of input files cannot be greater than 2")
    elif len(options.input) == 0:
        raise RuntimeError("BUG! Empty options.input array")

    if options.se:
        raise NotImplementedError("Sorry. Single-end alignments not implemented")

    return options

def main(argv=None):
    options = parse_args(argv)
    plugin = pyrapi.load_aligner('rapi_bwa')

    opts = plugin.opts()
    plugin.init(opts)
    _log.info("Using the %s aligner plugin, aligner version %s, plugin version %s",
            plugin.aligner_name(), plugin.aligner_version(), plugin.plugin_version())
    pe = not options.se
    _log.info("Working in %s mode", 'paired-end' if pe else 'single-end')
    _log.info("Reading data in %s format", options.format)

    _log.info("Loading reference %s", options.ref)
    ref = plugin.ref(options.ref)
    _log.info("Reference loaded")

    batch = plugin.read_batch(2 if pe else 1)

    input_generator = create_input(options)

    _log.info('loading data')
    for reads in input_generator:
        read1 = reads[0]
        batch.append(read1['id'], read1['seq'], read1['q'], plugin.QENC_SANGER)
        if pe:
            read2 = reads[1]
            batch.append(read2['id'], read2['seq'], read2['q'], plugin.QENC_SANGER)

    _log.info('Loaded %s reads', len(batch))

    _log.info("aligning...")
    aligner = plugin.aligner(opts)
    aligner.align_reads(ref, batch)
    _log.info("finished aligning")

    _log.info("Here's the output")
    for idx, fragment in enumerate(batch):
        _log.info('fragment')
        assert len(fragment) == 2
        sam = plugin.format_sam(batch, idx)
        print sam

    ref.unload()

if __name__ == '__main__':
    main()

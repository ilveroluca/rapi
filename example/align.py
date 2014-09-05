#!/usr/bin/env python

import argparse
import fileinput
import logging

import pyrapi.rapi as rapi

logging.basicConfig(level=logging.INFO)
_log = logging.getLogger('align')

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('ref')
    parser.add_argument('input', nargs='*')
    parser.add_argument('--se', action="store_true", default=False)
    options = parser.parse_args(args)
    return options

def main(argv=None):
    options = parse_args(argv)
    rapi.init(rapi.opts())

    _log.info("Loading reference %s", options.ref)
    ref = rapi.ref(options.ref)
    _log.info("Reference loaded")

    pe = not options.se
    _log.info("in %s mode", 'paired-end' if pe else 'single-end')

    batch = rapi.read_batch(2 if pe else 1)

    _log.info('input: %s', options.input)
    _log.info('reading data')
    for line in fileinput.input(files=options.input):
        row = line.rstrip('\n').split('\t')
        batch.append(row[0], row[1], row[2], rapi.QENC_SANGER)
        if pe:
            batch.append(row[0], row[3], row[4], rapi.QENC_SANGER)

    _log.info('Loaded %s reads', len(batch))

    ref.unload()

if __name__ == '__main__':
    main()

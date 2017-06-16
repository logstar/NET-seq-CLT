#!/usr/bin/env python2.7
import sys
import coordcount

def reverse_btcnt(ref_fn, ifn, ofn):
    bowtie_coord_count_table = coordcount.CoordCountTable(ifn, ref_fn)
    bowtie_coord_count_table.to_strand_reverse_count_tbl(ofn)

def main():
    argv = sys.argv
    if len(argv) != 4:
        sys.stderr.write('Usage:\n%s <reference fasta file> <input Bowtie count file (0-based)> \
<output strand reversed count file (0-based)>\n' % argv[0])
        return -2

    ref_fn = argv[1]
    ifn = argv[2]
    ofn = argv[3]

    reverse_btcnt(ref_fn, ifn, ofn)


if __name__ == "__main__":
    main()
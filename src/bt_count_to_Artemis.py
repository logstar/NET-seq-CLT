#!/usr/bin/env python2.7
import sys
import coordcount

def btcnt_to_artemis_cnt(ref_fn, ifn, ofn):
    bowtie_coord_count_table = coordcount.CoordCountTable(ifn, ref_fn)
    bowtie_coord_count_table.to_artemis_count_tbl(ofn)

def main():
    argv = sys.argv
    if len(argv) != 4:
        sys.stderr.write('Usage:\n%s <reference fasta file> <input Bowtie count file (0-based)> \
<output Artemis count file (1-based)>\n' % argv[0])
        return -2

    ref_fn = argv[1]
    ifn = argv[2]
    ofn = argv[3]
    btcnt_to_artemis_cnt(ref_fn, ifn, ofn)

if __name__ == "__main__":
    main()
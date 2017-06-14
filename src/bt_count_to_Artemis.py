#!/usr/bin/env python2.7
import sys
import coordcount

def btcnt_to_artemis_cnt(ref_fn, ifn, reverse_strand, ofn):
    bowtie_coord_count_table = coordcount.CoordCountTable(ifn, ref_fn)
    bowtie_coord_count_table.to_artemis_count_tbl(ofn, reverse_strand)

def main():
    argv = sys.argv
    if len(argv) != 5:
        sys.stderr.write('Usage:\n%s <reference fasta file> <input Bowtie count file (0-based)> \
<reverse strand?> <output Artemis count file (1-based)>\n' % argv[0])
        return -2

    ref_fn = argv[1]
    ifn = argv[2]
    reverse_strand = argv[3]

    if reverse_strand == 'reverse':
        reverse_strand = True
    elif reverse_strand == 'no-reverse':
        reverse_strand = False
    else:
        sys.stderr.write("<reverse strand?> can only be reverse or no-reverse\n")
        return -1
    
    ofn = argv[4]
    btcnt_to_artemis_cnt(ref_fn, ifn, reverse_strand, ofn)

if __name__ == "__main__":
    main()
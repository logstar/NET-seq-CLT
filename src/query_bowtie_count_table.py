#!/usr/bin/env python2.7
import sys
import coordcount

NUM_UP_BP = 50
NUM_DOWN_BP = 50

def is_valid_strand(strand):
    if strand in ('+', '-'):
        return True
    else:
        return False

def parse_coord_strand_file(ifn):
    coord_strand_list = []
    with open(ifn) as ifile:
        for line in ifile:
            fields = line.strip().split()
            coord = int(fields[0])

            if is_valid_strand(fields[1]):
                strand = fields[1]
            else:
                raise ValueError("<strand> can only be + or -")
            coord_strand_list.append((coord, strand))

    return tuple(coord_strand_list)

def main():
    argv = sys.argv
    if len(argv) != 5:
        sys.stderr.write('Usage:\n%s <reference fasta file> <input Bowtie count file (0-based)> \
<coordinate (0-based) strand list> <output seq count file (0-based)>\n' % argv[0])
        return -2

    ref_fn = argv[1]
    ifn = argv[2]
    cs_fn = argv[3]
    ofn = argv[4]

    coord_strand_list = parse_coord_strand_file(cs_fn)
    bowtie_coord_count_table = coordcount.CoordCountTable(ifn, ref_fn)

    ofile = open(ofn, 'w')
    ofile.write('\t'.join(['coord', 'strand', 'seq'] + ['rc_U' + str(i) for i in xrange(NUM_UP_BP, 0, -1)] + ['rc_0'] + ['rc_D' + str(i) for i in xrange(1, NUM_DOWN_BP + 1)]) + '\n')
    for cs_tup in coord_strand_list:
        coord = cs_tup[0]
        strand = cs_tup[1]
        seq_cnt_tup = bowtie_coord_count_table.get_up_down_stream_seq_cnt_tuple(coord, strand, NUM_UP_BP, NUM_DOWN_BP)
        ofile.write('\t'.join((str(coord), strand, seq_cnt_tup[0])) + '\t' + '\t'.join(map(str, seq_cnt_tup[1])) + '\n')

    ofile.close()

if __name__ == "__main__":
    main()

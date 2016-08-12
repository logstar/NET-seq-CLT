#!/usr/bin/env python2.7
import sys
import coordcount as cc
import collections
import fastautil

TssRecord = collections.namedtuple('TssRecord', 
                                   ['prmt_name', 'b_num', 'seq', 
                                    'strand', 'tss_coord'])
# Input coordinate is per strand 5' -> 3' coordinate 1-based.
# Parsed coordinate is forward strand 5' -> 3' coordinate 0-based.
def parse_tss_coord_strand_table(in_tss_fn, ref_fn):
    ref_genome = fastautil.SingleEntryFastaFile(ref_fn)

    tss_rec_list = []
    in_tss_file = open(in_tss_fn, 'r')
    for line in in_tss_file:
        fields = line.strip('\n').split('\t')
        prmt_name = fields[0]
        b_num = fields[1]
        seq = fields[2]
        if fields[3] == 'Forward':
            strand = '+'
            tss_coord = int(fields[4]) - 1
        elif fields[3] == 'Reverse':
            strand = '-'
            tss_coord = ref_genome.ref_seq_len - int(fields[4])
        else:
            raise ValueError('unknown strand %s' % fields[3])

        tss_rec_list.append(TssRecord(prmt_name, b_num, seq, strand, tss_coord))

    in_tss_file.close()
    return tss_rec_list

def get_tss_count_table(num_up_bp, num_down_bp, in_coord_cnt_fn, in_ref_fn, 
                        in_tss_fn, out_tss_cnt_fn):
    tss_rec_list = parse_tss_coord_strand_table(in_tss_fn, in_ref_fn)

    coord_count_table = cc.CoordCountTable(in_coord_cnt_fn, in_ref_fn)

    out_tss_cnt_file = open(out_tss_cnt_fn, 'w')
    for tss_rec in tss_rec_list:
        seq, count_tuple = coord_count_table.get_up_down_stream_seq_cnt_tuple(
            tss_rec.tss_coord, tss_rec.strand, num_up_bp, num_down_bp)

        min_seq_len = min(num_down_bp + 1, len(tss_rec.seq))
        assert seq[num_up_bp:num_up_bp + min_seq_len].upper() == tss_rec.seq[:min_seq_len]

        count_str_tuple = tuple(map(str, count_tuple))
        out_tss_cnt_file.write('\t'.join((tss_rec.prmt_name, tss_rec.b_num, 
            seq, tss_rec.strand, str(tss_rec.tss_coord)) + count_str_tuple) + '\n')
        
    out_tss_cnt_file.close()


def main():
    argv = sys.argv
    if len(argv) != 7:
        sys.stderr.write('Usage:\n%s <number of upstream bp> \
<number of downstream bp> <input coord count table> \
<input reference genome> <input TSS coord strand table> \
<output TSS count table>\n' % argv[0])
        return -2

    num_up_bp = int(argv[1])
    num_down_bp = int(argv[2])
    
    in_coord_cnt_fn = argv[3]
    in_ref_fn = argv[4]
    in_tss_fn = argv[5]
    out_tss_cnt_fn = argv[6]

    get_tss_count_table(num_up_bp, num_down_bp, in_coord_cnt_fn, in_ref_fn, 
                        in_tss_fn, out_tss_cnt_fn)

if __name__ == '__main__':
    main()
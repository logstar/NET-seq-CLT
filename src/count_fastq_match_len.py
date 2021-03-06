#!/usr/bin/env python2.7
import argparse
import dnautil
import sequtil
import fastqutil
import collections

# original seq is used to keep the original sequence without reverse / complement. 
TargetSeqRecord = collections.namedtuple('FastaQuery', ['coord', 'num_up_bp', 'num_down_bp', 
                                                        'strand', 'seq', 'original_seq'])

# reverse and complement must be boolean
def parse_fa_query_result_file(fa_qres_fn, reverse, complement):
    assert isinstance(reverse, bool)
    assert isinstance(complement, bool)

    tseq_rec_list = []
    with open(fa_qres_fn, 'r') as fa_qres_file:
        for line in fa_qres_file:
            fields = line.strip().split()
            assert len(fields) == 5
            coord = int(fields[0])
            num_up_bp = int(fields[1])
            num_down_bp = int(fields[2])
            strand = fields[3]
            original_seq = fields[4]
            target_seq = fields[4]

            if reverse:
                target_seq = target_seq[::-1]

            if complement:
                target_seq = dnautil.complement(target_seq)

            tseq_rec = TargetSeqRecord(coord = coord, num_up_bp = num_up_bp, 
                                      num_down_bp = num_down_bp, 
                                      strand = strand, seq = target_seq, 
                                      original_seq = original_seq)

            tseq_rec_list.append(tseq_rec)

    return tseq_rec_list


def count_fastq_match_len(fq_fn, target_seq_set):
    assert len(target_seq_set) == len(set(target_seq_set))
    
    tseq_mlen_cnter_list = map(lambda tseq: sequtil.TargetSeqMatchLenCounter(tseq), target_seq_set)

    for seqid, seq, rseqid, qscore in fastqutil.iterate_fastq_file(fq_fn):
        for tseq_mlen_cnter in tseq_mlen_cnter_list:
            tseq_mlen_cnter.add_query_seq(seq)

    return tseq_mlen_cnter_list

def all_same_mlen_range(tseq_mlen_cnt_list):
    if len(tseq_mlen_cnt_list) == 0:
        return True

    mlen_tup = tuple(map(lambda t: t[0], tseq_mlen_cnt_list[0][1:]))

    # assert all matched length ranges are the same
    for tseq_mlen_cnt in tseq_mlen_cnt_list:
        tseq_mlen_tup = tuple(map(lambda t: t[0], tseq_mlen_cnt[1:]))
        if tseq_mlen_tup != mlen_tup:
            return False

    return True

def combine_tseq_rec_tseq_mlen_cnt(tseq_rec, tseq_mlen_cnt):
    coord = tseq_rec.coord
    strand = tseq_rec.strand
    seq = tseq_rec.seq
    tseq = tseq_mlen_cnt[0]

    assert seq.upper() == tseq.upper()

    out_seq = tseq_rec.original_seq

    mlen_cnt_tup = tuple(map(lambda t: t[1], tseq_mlen_cnt[1:]))
    
    return (coord, strand, out_seq) + mlen_cnt_tup

# tseq_rec_list is a list of TargetSeq
# tseq_mlen_cnt_list is a list of tseq_mlen_cnt = [seq, (mlen0, mlen0_cnt), 
#                                                  (mlen1, mlen1_cnt), (mlen0, mlen2_cnt), ...]
# reverse or complement output sequence
def fmt_tseq_rec_mlen_cnt(tseq_rec_list, tseq_mlen_cnt_list):
    if len(tseq_rec_list) == 0:
        raise ValueError('ftseq_rec_list is empty.')
    
    if len(tseq_mlen_cnt_list) == 0:
        raise ValueError('tseq_mlen_cnt_list is empty.')
    
    assert all_same_mlen_range(tseq_mlen_cnt_list), 'matched length ranges are different'
    mlen_tup = tuple(map(lambda t: t[0], tseq_mlen_cnt_list[0][1:]))
    
    tseq_rec_mlen_cnt_list = []

    tseq_list = map(lambda l: l[0].upper(), tseq_mlen_cnt_list)

    assert len(tseq_list) == len(set(tseq_list)), 'target sequences are duplicated'

    for tseq_rec in tseq_rec_list:
        if tseq_rec.seq.upper() not in tseq_list:
            raise ValueError('Cannot find target seq.')

        for tseq_mlen_cnt in tseq_mlen_cnt_list:
            if tseq_rec.seq.upper() == tseq_mlen_cnt[0].upper():
                tseq_rec_mlen_cnt_list.append(
                    combine_tseq_rec_tseq_mlen_cnt(tseq_rec, tseq_mlen_cnt))

    assert len(tseq_rec_mlen_cnt_list) > 0
    ncol = map(lambda t: len(t), tseq_rec_mlen_cnt_list)
    assert len(set(ncol)) == 1
    assert ncol[0] == len(mlen_tup) + 3
    
    header = 'coord\tstrand\tseq\t' + '\t'.join(map(str, mlen_tup))
    rows = map(lambda tup: '\t'.join(map(str, tup)), tseq_rec_mlen_cnt_list)

    tbl_str = header + '\n' + '\n'.join(rows) + '\n'

    return tbl_str


def main():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('fa_qres_fn', metavar = '<fasta query result file>',
                            help = 'A tab delimited file listing all fasta query results.' 
                                   'These are target sequences for matching FASTQ reads.'
                                   'Fields from left to right are:'
                                   'coordinate, number of upstream bp, '
                                   'number of downstream bp, strand, and sequence')

    arg_parser.add_argument('-r', '--reverse', 
                            help = 'Sepcify this argument to reverse the FASTA'
                                   'target sequences before calculating matched length.'
                                   'NOTE: this argument does not affect output sequence.',
                            action = 'store_true')

    arg_parser.add_argument('-c', '--complement', 
                            help = 'Sepcify this argument to complement the FASTA'
                                   'target sequences before calculating matched length.'
                                   'NOTE: this argument does not affect output sequence.',
                            action = 'store_true')

    arg_parser.add_argument('fq_fn', metavar = '<fastq file for finding the max len>')

    arg_parser.add_argument('o_fn', metavar = '<output file>')
    
    args = arg_parser.parse_args()
    
    tseq_rec_list = parse_fa_query_result_file(args.fa_qres_fn, reverse = args.reverse, 
                                               complement = args.complement)

    target_seq_set = set(map(lambda tseq_rec: tseq_rec.seq.upper(), tseq_rec_list))

    tseq_mlen_cnter_list = count_fastq_match_len(args.fq_fn, target_seq_set)

    tseq_mlen_cnt_list = map(lambda cnter: cnter.to_list(0, 5), tseq_mlen_cnter_list)

    o_tbl_str = fmt_tseq_rec_mlen_cnt(tseq_rec_list, tseq_mlen_cnt_list)

    with open(args.o_fn, 'w') as ofile:
        ofile.write(o_tbl_str)


if __name__ == '__main__':
    main()
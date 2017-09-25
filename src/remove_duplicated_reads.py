#!/usr/bin/env python2.7
from __future__ import print_function
import argparse
import fastqutil

def remove_duplicated_reads(num_qc_bp, qcutoff, dc_len, ifn, ofn):
    dc_seq_set = set()
    short_read_cnter = {}
    
    o_fastq_file = open(ofn, 'w')
    for seqid, seq, rseqid, qscore in fastqutil.iterate_fastq_file(ifn):
        # If the read is shorter than dc_len (dup checking len), discard the read.
        if len(seq) < dc_len:
            short_read_cnter[len(seq)] = short_read_cnter.get(len(seq), 0) + 1
        else:
            dc_seq = seq[:dc_len]
            if ((min(qscore[:num_qc_bp]) >= qcutoff) and ('N' not in seq[:num_qc_bp]) 
                and (dc_seq not in dc_seq_set)):
                dc_seq_set.add(dc_seq)
                o_fastq_file.write('\n'.join((seqid, seq, rseqid, qscore)) + '\n')

    o_fastq_file.close()

    print('Sequencing reads shorther than {}'.format(dc_len))
    print('Length\tCount')
    for short_read_len in sorted(short_read_cnter.keys()):
        print('%d\t%d'.format(short_read_len, short_read_cnter[short_read_cnter]))
        
    return

def main():
    arg_parser = argparse.ArgumentParser()

    arg_parser.add_argument('num_qc_bp', type = int,
                            metavar = '<number of quality checking bases>')

    arg_parser.add_argument('qc_int', type = int,
                            metavar = '<Phred quality score cutoff>',
                            help = 'FASTQ per base Phred quality score cutoff.'
                                   'Discard reads with quality score < qcutoff')

    arg_parser.add_argument('dc_len', type = int,
                            metavar = '<length of sequence for duplicate checking>')

    arg_parser.add_argument('ifn', metavar = '<input fastq file>')

    arg_parser.add_argument('ofn', metavar = '<output fastq file>',
                            help = 'The first occurence of the read passing '
                                   'quality check.')

    args = arg_parser.parse_args()
    
    remove_duplicated_reads(args.num_qc_bp, fastqutil.int_to_qscore(args.qc_int), 
                            args.dc_len, args.ifn, args.ofn)

if __name__ == "__main__":
    main()

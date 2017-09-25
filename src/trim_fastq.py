#!/usr/bin/env python2.7
import argparse
import fastqutil

def trim_fastq(i_fn_list, trim_beg_ind, final_length, qcutoff, o_stat_fn, o_fastq_fn):
    length_failed_reads = 0
    qual_failed_reads = 0
    total_reads = 0
    filtered_reads = 0

    o_fastq_file = open(o_fastq_fn, 'w')
    for seq_fn in i_fn_list:
        for seqid, seq, rseqid, qscore in fastqutil.iterate_fastq_file(seq_fn):
            total_reads += 1
            if len(seq) >= trim_beg_ind + final_length:
                seq_trim = seq[trim_beg_ind : trim_beg_ind + final_length]
                qscore_trim = qscore[trim_beg_ind : trim_beg_ind + final_length]
                if min(qscore_trim) >= qcutoff:
                    filtered_reads += 1
                    o_fastq_file.write('\n'.join((seqid, seq_trim, rseqid, qscore_trim)) + '\n')
                else:
                    qual_failed_reads += 1
            else:
                length_failed_reads += 1
    o_fastq_file.close()

    stats = '\n'.join(i_fn_list) + '\n'
    stats += 'Total reads: %d\n' % total_reads
    stats += 'Quality score failed reads: %d\n' % qual_failed_reads
    stats += 'Length failed reads: %d\n' % length_failed_reads
    stats += 'Filtered reads: %d (%s%%)\n' % (filtered_reads, 
        "{:.4f}".format(float(filtered_reads) / total_reads * 100))

    with open(o_stat_fn, 'w') as o_stat_file:
        o_stat_file.write(stats)

    return

def main():
    arg_parser = argparse.ArgumentParser()

    arg_parser.add_argument('qc_int', type = int,
                            metavar = '<Phred quality score cutoff>',
                            help = 'FASTQ per base Phred quality score cutoff.'
                                   'Discard reads with quality score < qcutoff')

    arg_parser.add_argument('trim_beg_ind', type = int,
                            metavar = '<5\' trim length>')

    arg_parser.add_argument('final_length', type = int,
                            metavar = '<final sequence length>')

    arg_parser.add_argument('o_stat_fn', metavar = '<output stats>')

    arg_parser.add_argument('o_fastq_fn', metavar = '<output trimmed FASTQ file>')

    arg_parser.add_argument('i_fn_list', nargs = '+', 
                            metavar = '<input fastq files (separated by space)>')
    
    args = arg_parser.parse_args()

    trim_fastq(args.i_fn_list, args.trim_beg_ind, args.final_length, 
               fastqutil.int_to_qscore(args.qc_int), args.o_stat_fn, 
               args.o_fastq_fn)

if __name__ == "__main__":
    main()

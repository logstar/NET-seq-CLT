#!/usr/bin/env python2.7
import sys
import fastqutil

def remove_duplicated_reads(num_qc_bp, qcutoff, dc_len, ifn, ofn):
    dc_seq_set = set()

    o_fastq_file = open(ofn, 'w')
    for seqid, seq, rseqid, qscore in fastqutil.iterate_fastq_file(ifn):
        dc_seq = seq[:dc_len]
        if ((min(qscore[:num_qc_bp]) >= qcutoff) and ('N' not in seq[:num_qc_bp]) 
            and (dc_seq not in dc_seq_set)):
            dc_seq_set.add(dc_seq)
            o_fastq_file.write('\n'.join((seqid, seq, rseqid, qscore)) + '\n')
    o_fastq_file.close()
    return

def main():
    argv = sys.argv
    if len(argv) != 6:
        sys.stderr.write('Usage:\n%s <number of quality checking bases> <Phred quality score (>=)> \
<length of sequence for duplicate checking> <input fastq file> <output fastq file>\n\
The first occurence of the read passing quality check.\n' % argv[0])
        return -2

    num_qc_bp = int(argv[1])

    qc = int(argv[2])
    if qc < 0 or qc > 93:
        sys.stderr.write('Quality score cutoff should >= 0 and <= 93')
        return -1
    qcutoff = chr(qc + 33)

    dc_len = int(argv[3])
    ifn = argv[4]
    ofn = argv[5]

    remove_duplicated_reads(num_qc_bp, qcutoff, dc_len, ifn, ofn)

if __name__ == "__main__":
    main()

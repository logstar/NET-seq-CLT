#!/ingens/home/yuanchao/tools/miniconda2/envs/bowtie_alignment/bin/python2.7
import sys
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
    argv = sys.argv
    if len(argv) <= 6:
        sys.stderr.write('Usage:\n%s [Phred quality score cutoff (>=)] [5\' trim length] \
[final sequence length] <output stats> <output trimmed FASTQ file> \
<input fastq files (separated by space)>\n' % argv[0])
        return -2

    qc = int(argv[1])
    if qc < 0 or qc > 93:
        sys.stderr.write('Quality score cutoff should >= 0 and <= 93')
        return -1

    qcutoff = chr(qc + 33)

    trim_beg_ind = int(argv[2]) # 0-based index, so trim length = index of the first trimmed base
    final_length = int(argv[3])

    o_stat_fn = argv[4]
    o_fastq_fn = argv[5]

    i_fn_list = argv[6:]

    trim_fastq(i_fn_list, trim_beg_ind, final_length, qcutoff, o_stat_fn, o_fastq_fn)

if __name__ == "__main__":
    main()

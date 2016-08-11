#!/usr/bin/env python2.7
import sys
import fastautil
import bowtieutil as bu

def parse_bt_file(in_bt_fn, qcutoff, ref_length):
    num_uniq_aligned = 0
    num_qual_failed = 0

    bt_align_counter = bu.BowtieRecordCounter(ref_length)
    #{position : {'+' : count, '-' : count}, ...}
    for bt_rec_ntup in bu.iterate_bowtie_out_file(in_bt_fn):
        num_uniq_aligned += 1
        assert bt_rec_ntup.mm_desc == '', ("Read alignment has mismatches. %s" % 
                                           bt_rec_ntup._asdict())

        if min(bt_rec_ntup.qscore) < qcutoff:
            num_qual_failed += 1
        else:
            bt_align_counter.insert_bt_rec_ntup(bt_rec_ntup)
    return [bt_align_counter, num_uniq_aligned, num_qual_failed]

def write_parse_bt_stats(stats_list, stats_fn):
    num_uniq_aligned = stats_list[0]
    num_qual_failed = stats_list[1]
    num_filtered = num_uniq_aligned - num_qual_failed

    stats_file = open(stats_fn, 'w')
    stats_file.write("Number of reads uniquely aligned:%d\n" % num_uniq_aligned)
    stats_file.write("Number of reads with quality socre <= cutoff: %d (%s%%)\n" 
                      % (num_qual_failed, "{:.4f}".format(float(num_qual_failed) / num_uniq_aligned * 100)))
    stats_file.write("Number of uniquely aligned reads passed quality socre filter: %d (%s%%)\n" 
                      % (num_filtered, "{:.4f}".format(float(num_filtered) / num_uniq_aligned * 100)))

    stats_file.close()

def main():
    argv = sys.argv
    if len(argv) != 6:
        sys.stderr.write('Usage:\n%s [>= quality score] <input reference genome file name> \
<input bowtie output> <output stats file> <output start position count file 0 based>\n' % argv[0])
        return -2

    qc = int(argv[1])
    if qc < 0 or qc > 93:
        sys.stderr.write('Quality score cutoff should >= 0 and <= 93')
        return -1

    qcutoff = chr(qc + 33)

    in_ref_fn = argv[2]
    in_ref_fa_file = fastautil.SingleEntryFastaFile(in_ref_fn)
    ref_length = in_ref_fa_file.ref_seq_len

    in_bt_fn = argv[3]
    result_list = parse_bt_file(in_bt_fn, qcutoff, ref_length)

    ostats_fn = argv[4]
    write_parse_bt_stats(result_list[1:], ostats_fn)

    ocount_fn = argv[5]
    result_list[0].output_count_table(ocount_fn)


if __name__ == "__main__":
    main()

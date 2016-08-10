#!/usr/bin/env python2.7
import sys
def parse_bt_file(in_bt_fn, qcutoff, ref_length):
#    ref_length = 4641652
    trim_length = 30
    num_uniq_aligned = 0
    num_qual_failed = 0

    align_dict = {}
    #{position : [+ count, - count]}
    in_bt_file = open(in_bt_fn, 'r')
    for line in in_bt_file:
        num_uniq_aligned += 1

        fields = line[:-1].split('\t')
        mismatch_des = fields[-1]
        num_other_align = fields[-2]
        if mismatch_des != '' or num_other_align != '0':
            print line
            sys.exit(2)

        qual_score = fields[5]
        qual_failed = False
        for q in qual_score:
            if q < qcutoff:
                num_qual_failed += 1
                qual_failed = True
                break
        if qual_failed:
            continue

        strand = fields[1]
        if strand == '+':
            position = int(fields[3])
            cl_ind = 0
        else:
            position = int(fields[3]) + trim_length - 1
            if position >= ref_length:
                position = position - ref_length
                print 1
            cl_ind = 1
        
        if position in align_dict:
            align_dict[position][cl_ind] = align_dict[position][cl_ind] + 1
        else:
            align_dict[position] = [0, 0]
            align_dict[position][cl_ind] = 1
    
    in_bt_file.close()
    return [align_dict, num_uniq_aligned, num_qual_failed]

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

def write_align_dict(align_dict, ocount_fn):
    ocount_file = open(ocount_fn, 'w')

    for key in sorted(align_dict.keys()):
        ocount_file.write("%d\t%d\t%d\n" % (key, align_dict[key][0], align_dict[key][1]))

    ocount_file.close()


def get_ref_genome_length(ref_genome_fn):
    ref_genome_file = open(ref_genome_fn)

    description = ref_genome_file.readline()[:-1]
    ref_seq_line = ref_genome_file.readline()[:-1]
    ref_seq = ''
    while ref_seq_line != '':
        ref_seq = ref_seq + ref_seq_line
        ref_seq_line = ref_genome_file.readline()[:-1]

    ref_genome_file.close()

    ref_len = len(ref_seq)

    return ref_len



if __name__ == "__main__":
    argv = sys.argv
    if len(argv) != 6:
        sys.stderr.write('Usage:\n%s [>= quality score] <input reference genome file name> \
<input bowtie output> <output stats file> <output start position count file 0 based>\n' % argv[0])
        sys.exit(2)

    qc = int(argv[1])
    if qc < 0 or qc > 93:
        sys.stderr.write('Quality score cutoff should >= 0 and <= 93')
        sys.exit(2)

    asciiseq = "!\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
    qcutoff = asciiseq[qc]

    in_ref_fn = argv[2]
    ref_length = get_ref_genome_length(in_ref_fn)

    in_bt_fn = argv[3]
    result_list = parse_bt_file(in_bt_fn, qcutoff, ref_length)

    ostats_fn = argv[4]
    write_parse_bt_stats(result_list[1:], ostats_fn)

    ocount_fn = argv[5]
    write_align_dict(result_list[0], ocount_fn)

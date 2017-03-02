#!/usr/bin/env python2.7
import sys

def get_tss_read_count_dict(in_fn):
	tss_read_count_dict = {}
	with open(in_fn, 'r') as in_file:
		for line in in_file:
			fields = line[:-1].split()
			position_b0 = int(fields[0])
			pos_count = int(fields[1])
			neg_count = int(fields[2])

			if position_b0 in tss_read_count_dict:
				sys.stderr.write("Position duplicated: %d\n" % position_b0)
			else:
				tss_read_count_dict[position_b0] = [pos_count, neg_count]
	return tss_read_count_dict

def merge_tss_read_count_dicts(tss_read_count_dict_list):
	merged_trc_dict = {}
	for trc_dict in tss_read_count_dict_list:
		for key in trc_dict:
			if key not in merged_trc_dict:
				merged_trc_dict[key] = [0, 0]

			merged_trc_dict[key][0] += trc_dict[key][0]
			merged_trc_dict[key][1] += trc_dict[key][1]

	return merged_trc_dict

def write_align_dict(align_dict, ocount_fn):
    ocount_file = open(ocount_fn, 'w')

    for key in sorted(align_dict.keys()):
        ocount_file.write("%d\t%d\t%d\n" % (key, align_dict[key][0], align_dict[key][1]))

    ocount_file.close()

if __name__ == "__main__":
	argv = sys.argv
	if len(argv) <= 3:
		sys.stderr.write("Usage:\n%s <output merged TSS read count file>\
<input TSS read count files (0 based)>\nInput files separated by space.\n" % argv[0])
		sys.exit(2)

	o_fn = argv[1]
	in_fn_list = argv[2:]
	tss_read_count_dict_list = map(get_tss_read_count_dict, in_fn_list)
	print len(tss_read_count_dict_list)
	merged_tss_read_count_dict = merge_tss_read_count_dicts(tss_read_count_dict_list)
	write_align_dict(merged_tss_read_count_dict, o_fn)



#!/usr/bin/env python2.7
import sys
from merge_tss_read_count_files import get_tss_read_count_dict

def get_tss_ref_list(in_tss_ref_fn):
	tss_ref_list = []
	with open(in_tss_ref_fn, 'r') as in_tss_ref_file:
		for line in in_tss_ref_file:
			fields = line[:-1].split('\t')
			coord = int(fields[0])
			strand = fields[1]
			seq = fields[2]
			nearest_genes = fields[-4:]

			tss_ref_list.append(tuple([coord, strand, seq] + nearest_genes))

	return tss_ref_list

def gen_pct_tss_file(in_tss_ct_fn, rc_cutoff, tss_ref_list, o_stats_fn, o_tss_pct_fn):
	ref_length = 4641652

	tss_rc_dict = get_tss_read_count_dict(in_tss_ct_fn)

	pos_rc_list = [0] * ref_length
	neg_rc_list = [0] * ref_length

	for tss_coord in tss_rc_dict:
		if pos_rc_list[tss_coord] != 0 or neg_rc_list[tss_coord] != 0:
			raise ValueError("Duplicated TSS coordinate.")
		pos_rc_list[tss_coord] = tss_rc_dict[tss_coord][0]
		neg_rc_list[tss_coord] = tss_rc_dict[tss_coord][1]

	tss_pct_list = []

	num_pos_ref_tss = 0
	num_neg_ref_tss = 0
	num_pos_arcc_tss = 0
	num_neg_arcc_tss = 0
	
	for tss_ref_tup in tss_ref_list:
		tss_coord_b0 = tss_ref_tup[0]
		strand = tss_ref_tup[1]
		window_range = range(tss_coord_b0 - 5, tss_coord_b0 + 6)

		if tss_coord_b0 + 5 >= ref_length:
			for i in xrange(len(window_range)):
				if window_range[i] >= ref_length:
					window_range[i] = window_range[i] - ref_length

		if strand == '+':
			num_pos_ref_tss += 1
			
			pos_count_window = [pos_rc_list[i] for i in window_range]
			pos_count_window_sum = sum(pos_count_window)

			if pos_count_window_sum >= rc_cutoff:
				num_pos_arcc_tss += 1
				pos_pct_window = [str(pos_rc / float(pos_count_window_sum) * 100) 
								  for pos_rc in pos_count_window]
			else:
				pos_pct_window = [''] * 11

			tss_pct_list.append(tuple([str(tss_coord_b0)] + list(tss_ref_tup[1:3]) + 
									   map(str, pos_count_window) + [str(pos_count_window_sum)] + 
									   pos_pct_window + list(tss_ref_tup[3:])))

		if strand == '-':
			num_neg_ref_tss += 1

			#Ref coordinates reversed, so that list ordered from upstream to downstream
			neg_count_window = [neg_rc_list[i] for i in window_range[::-1]]
			neg_count_window_sum = sum(neg_count_window)

			if neg_count_window_sum >= rc_cutoff:
				num_neg_arcc_tss += 1
				neg_pct_window = [str(neg_rc / float(neg_count_window_sum) * 100) 
								  for neg_rc in neg_count_window]
			else:
				neg_pct_window = [''] * 11

			tss_pct_list.append(tuple([str(tss_coord_b0)] + list(tss_ref_tup[1:3]) + 
									   map(str, neg_count_window) + [str(neg_count_window_sum)] + 
									   neg_pct_window + list(tss_ref_tup[3:])))


	with open(o_tss_pct_fn, 'w') as o_tss_pct_file:
		for tss_pct_tup in tss_pct_list:
			o_tss_pct_file.write('\t'.join(tss_pct_tup) + '\n')

	stats = "# reference TSSs on + strand: %d \n" % num_pos_ref_tss
	stats += "# TSSs on + strand passed read count filter: %d (%s%%) \n" % (num_pos_arcc_tss,
				"{:.4f}".format(float(num_pos_arcc_tss) / num_pos_ref_tss * 100))
	stats += "# reference TSSs on - strand: %d \n" % num_neg_ref_tss
	stats += "# TSSs on - strand passed read count filter: %d (%s%%) \n" % (num_neg_arcc_tss,
				"{:.4f}".format(float(num_neg_arcc_tss) / num_neg_ref_tss * 100))

	with open(o_stats_fn, 'w') as o_stats_file:
		o_stats_file.write(stats)




if __name__ == "__main__":
	argv = sys.argv
	if len(argv) != 6:
		sys.stderr.write("Usage:\n%s <input TSS reference file (0 based)>\
 <input TSS count file> [read count cutoff (>=)] <output stats file>\
 <output percent TSS file> \n" % argv[0])
		sys.exit(2)
	in_tss_ref_fn = argv[1]
	in_tss_ct_fn = argv[2]
	rc_cutoff = int(argv[3])
	o_stats_fn = argv[4]
	o_tss_pct_fn = argv[5]

	tss_ref_list = get_tss_ref_list(in_tss_ref_fn)

	gen_pct_tss_file(in_tss_ct_fn, rc_cutoff, tss_ref_list, o_stats_fn, o_tss_pct_fn)

#!/usr/bin/env python2.7
import sys
RC_CUTOFF = 20

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
				sys.exit(2)
			else:
				tss_read_count_dict[position_b0] = [pos_count, neg_count]
	return tss_read_count_dict

def rev_comp(seq):
        refd = dict(zip('ACGTacgt','TGCAtgca'))
        cseq = ''.join(map(lambda i: refd[i], seq))
        rcseq = cseq[::-1]
        return rcseq

def str_read_genome(ref_genome_file):
    description = ref_genome_file.readline()[:-1]
    ref_seq_line = ref_genome_file.readline()[:-1]
    ref_seq = ''
    while ref_seq_line != '':
        ref_seq = ref_seq + ref_seq_line
        ref_seq_line = ref_genome_file.readline()[:-1]

    return ref_seq

def get_tss_list(tss_read_count_dict, ref_str):
	num_pos_coords_have_count = 0
	sum_pos_rc = 0
	num_neg_coords_have_count = 0
	sum_neg_rc = 0
	ref_length = len(ref_str)
	pos_rc_list = [0] * ref_length
	neg_rc_list = [0] * ref_length

	for key in tss_read_count_dict:
		position_b0 = key
		pos_count = tss_read_count_dict[key][0]
		neg_count = tss_read_count_dict[key][1]
		if pos_count != 0:
			num_pos_coords_have_count += 1
			sum_pos_rc += pos_count
		if neg_count != 0:
			num_neg_coords_have_count += 1
			sum_neg_rc += neg_count

		pos_rc_list[position_b0] += pos_count
		neg_rc_list[position_b0] += neg_count

	tss_list = []

	position_b0 = 0
	pos_count = 0
	neg_count = 0

	num_pos_arc_coords = 0 #Number of positions on + strand with read count >= Read count cutoff
	num_pos_tie_coords = 0 
	#Number of Positions on + strand are local maxima and have local ties within 11bp
	num_pos_lm = 0 #Number of positions on + strand with RC >= 50 and no higher RC within 11bp
	num_pos_tss = 0 #Number of TSSs on + strand
	num_neg_arc_coords = 0
	num_neg_tie_coords = 0
	num_neg_lm = 0
	num_neg_tss = 0
	for position_b0 in xrange(ref_length):
		pos_count = pos_rc_list[position_b0]
		neg_count = neg_rc_list[position_b0]
		window_range = range(position_b0 - 5, position_b0) + range(position_b0 + 1, position_b0 + 6)
		if position_b0 + 5 >= ref_length:
			for i in xrange(len(window_range)):
				if window_range[i] >= ref_length:
					window_range[i] = window_range[i] - ref_length

		if pos_count >= RC_CUTOFF:
			num_pos_arc_coords += 1
			is_lm = True #is local maxima
			has_lt = False #has local tie
			for local_position_b0 in window_range:
				if pos_count < pos_rc_list[local_position_b0]:
					is_lm = False
					break
				if pos_count == pos_rc_list[local_position_b0]:
					has_lt = True

			if is_lm:
				num_pos_lm += 1
				is_tss = True #is TSS
				if has_lt:
					num_pos_tie_coords += 1
					for local_position_b0 in window_range[:5]:
						if pos_count == pos_rc_list[local_position_b0]:
							is_tss = False
				if is_tss:
					num_pos_tss += 1
					seq_start = position_b0 - 50
					seq_end = position_b0 + 50
					if seq_start < 0:
						seq_start = 0
					if seq_end >= ref_length:
						seq_end = ref_length - 1
					up_seq = ref_str[seq_start : position_b0] 
					tss_seq = ref_str[position_b0].lower() 
					down_seq = ref_str[position_b0 + 1 : seq_end + 1]
					seq = up_seq + tss_seq + down_seq
					pos_count_window_range = window_range[:5] + [position_b0] + window_range[5:]
					pos_count_window = [str(pos_rc_list[i]) for i in pos_count_window_range]

					tss_list.append('\t'.join([str(position_b0), '+', seq] + pos_count_window))
	
		if neg_count >= RC_CUTOFF:
			num_neg_arc_coords += 1
			is_lm = True
			has_lt = False
			for local_position_b0 in window_range:
				if neg_count < neg_rc_list[local_position_b0]:
					is_lm = False
					break
				if neg_count == neg_rc_list[local_position_b0]:
					has_lt = True

			if is_lm:
				num_neg_lm += 1
				is_tss = True
				if has_lt:
					num_neg_tie_coords += 1
					for local_position_b0 in window_range[:5]:
						if neg_count == neg_rc_list[local_position_b0]:
							is_tss = False
				if is_tss:
					num_neg_tss += 1
					seq_start = position_b0 - 50
					seq_end = position_b0 + 50
					if seq_start < 0:
						seq_start = 0
					if seq_end >= ref_length:
						seq_end = ref_length - 1
					up_seq = ref_str[seq_start : position_b0] 
					tss_seq = ref_str[position_b0].lower() 
					down_seq = ref_str[position_b0 + 1 : seq_end + 1]
					seq = up_seq + tss_seq + down_seq
					rc_seq = rev_comp(seq)
					neg_count_window_range = window_range[:5] + [position_b0] + window_range[5:]
					neg_count_window = [str(neg_rc_list[i]) for i in neg_count_window_range]
					rev_neg_count_window = neg_count_window[::-1]

					tss_list.append('\t'.join([str(position_b0), '-', rc_seq] + rev_neg_count_window))
		
	stats_str = "# Reads aligned to + strand: %d\n" % sum_pos_rc
	stats_str += "# coordinates in + strand have read count > 0: %d\n" % num_pos_coords_have_count
	stats_str += "# coordinates in + strand have read count >= %d: %d\n" % (RC_CUTOFF, num_pos_arc_coords)
	stats_str += "# local maxima in + strand (coordinates in + strand have read count >= %d \
with no higer read count in the local 11bp window): %d\n" % (RC_CUTOFF, num_pos_lm)
	stats_str += "# local maxima in + strand have local tie: %d\n" % num_pos_tie_coords
	stats_str += "# TSS in + strand: %d\n" % num_pos_tss

	stats_str += '--------------------------------------------------------\n'
	stats_str += "# Reads aligned to - strand: %d\n" % sum_neg_rc
	stats_str += "# coordinates in - strand have read count > 0: %d\n" % num_neg_coords_have_count
	stats_str += "# coordinates in - strand have read count >= %d: %d\n" % (RC_CUTOFF, num_neg_arc_coords)
	stats_str += "# local maxima in - strand (coordinates in - strand have read count >= %d \
with no higer read count in the local 11bp window): %d\n" % (RC_CUTOFF, num_neg_lm)
	stats_str += "# local maxima in - strand have local tie: %d\n" % num_neg_tie_coords
	stats_str += "# TSS in - strand: %d\n" % num_neg_tss

	return (tss_list, stats_str)






if __name__ == "__main__":
	argv = sys.argv
	if len(argv) != 5:
		sys.stderr.write("Usage:\n%s <input TSS read count file (0 based)>\
 <input reference genome FASTA> <output stats file> <output TSS file> \n" % argv[0])
		sys.exit(2)

	in_trc_fn = argv[1]
	in_ref_fn = argv[2]
	o_stats_fn = argv[3]
	o_tss_fn = argv[4]

	tss_read_count_dict = get_tss_read_count_dict(in_trc_fn)
	ref_str = str_read_genome(open(in_ref_fn, 'r'))
	result_tuple = get_tss_list(tss_read_count_dict, ref_str)

	with open(o_stats_fn, 'w') as o_stats_file:
		o_stats_file.write(result_tuple[1])

	with open(o_tss_fn, 'w') as o_tss_file:
		o_tss_file.write('\n'.join(result_tuple[0]))


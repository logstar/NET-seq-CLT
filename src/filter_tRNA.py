#!/usr/bin/env python2.7
import sys
from intervaltree import IntervalTree

def gen_tRNA_interval_tree(in_gff_fn, o_trna_gff_fn):
	trna_interval_tree = IntervalTree()

	in_gff_file = open(in_gff_fn, 'r')
	o_trna_gff_file = open(o_trna_gff_fn, 'w')
	print o_trna_gff_fn
	for line in in_gff_file:
		if line[0] == '#':
			continue
		fields = line[:-1].split('\t')
		if fields[2] == "tRNA":
			#seqname = fields[0]
			#source = fields[1]
			#feature = fields[2]
			start_b0 = int(fields[3]) - 1
			end_b0 = int(fields[4]) - 1
			#strand = fields[6]
			#attribute = fields[-1]

			trna_interval_tree[start_b0 : end_b0 + 0.5] = line
			o_trna_gff_file.write(line)

	in_gff_file.close()
	o_trna_gff_file.close()
	return trna_interval_tree

def filter_tRNA_reads(in_trc_fn, o_stats_fn, o_trc_fn, trna_interval_tree):
	num_total = 0
	num_tRNA = 0

	in_trc_file = open(in_trc_fn, 'r')
	o_trc_file = open(o_trc_fn, 'w')
	for line in in_trc_file:
		fields = line[:-1].split()
		position_b0 = int(fields[0])
		pos_count = int(fields[1])
		neg_count = int(fields[2])

		num_total = num_total + pos_count + neg_count

		overlap_trna_set = trna_interval_tree[position_b0]
		if len(overlap_trna_set) != 0:
			for trna_interal in overlap_trna_set:
				strand = trna_interal[2].split('\t')[6]
				if strand == '+':
					num_tRNA += pos_count
					pos_count = 0
				elif strand == '-':
					num_tRNA += neg_count
					neg_count = 0
				else:
					sys.stderr.write("Strand unexpected: %s" % strand)
					sys.exit(2)

		o_trc_file.write('\t'.join(map(str, [position_b0, pos_count, neg_count])) + '\n')

	o_trc_file.close()
	in_trc_file.close()

	num_filtered = num_total - num_tRNA
	o_stats_file = open(o_stats_fn, 'w')
	o_stats_file.write("Number of reads uniquely aligned:%d\n" % num_total)
	o_stats_file.write("Number of reads aligned to tRNA: %d (%s%%)\n" 
						% (num_tRNA, "{:.4f}".format(float(num_tRNA) / num_total * 100)))
	o_stats_file.write("Number of filtered reads: %d (%s%%)\n" 
						% (num_filtered, "{:.4f}".format(float(num_filtered) / num_total * 100)))

	o_stats_file.close()


if __name__ == "__main__":
	argv = sys.argv
	if len(argv) != 6:
		sys.stderr.write('Usage:\n%s <input GFF file 1 based> <input coord count file 0 based> \
<output filter stats> <output filtered TSS read count file 0 based> \
<output tRNA GFF>\n' % argv[0])
		sys.exit(2)

	in_gff_fn = argv[1]
	o_trna_gff_fn = argv[5]
	trna_interval_tree = gen_tRNA_interval_tree(in_gff_fn, o_trna_gff_fn)

	in_trc_fn = argv[2] #Tss Read Count file
	o_stats_fn = argv[3]
	o_trc_fn = argv[4]
	filter_tRNA_reads(in_trc_fn, o_stats_fn, o_trc_fn, trna_interval_tree)



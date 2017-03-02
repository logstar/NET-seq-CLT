from intervaltree import IntervalTree

class GFF(object):
	"""Simple GFF class for E coli genome. """
	
	def __init__(self, in_gff_fn):
		self.ref_length = -1
		self.trna_interval_tree = IntervalTree()
		self.gene_pos_interval_tree = IntervalTree()
		self.gene_neg_interval_tree = IntervalTree()
		self.gene_pos_beg_coord_dict = {}
		self.gene_pos_end_coord_dict = {}
		self.gene_neg_beg_coord_dict = {}
		self.gene_neg_end_coord_dict = {}

		in_gff_file = open(in_gff_fn, 'r')
		for line in in_gff_file:
			if line[0] == '#':
				if line[:17] == '##sequence-region':
					print line.strip().split(' ')
					self.ref_length = int(line.strip().split(' ')[-1])
					print "ref lengh = %d" % self.ref_length
				continue
			if self.ref_length == -1:
				in_gff_file.close()
				raise ValueError("Cannot find ##sequence-region line in gff file %s" % in_gff_fn)
			fields = line[:-1].split('\t')
			beg_b0 = int(fields[3]) - 1
			end_b0 = int(fields[4]) - 1
			strand = fields[6]

			if fields[2] == "tRNA":
				self.trna_interval_tree[beg_b0 : end_b0 + 0.5] = fields
			if fields[2] == "gene":
				if strand == '+':
					self.gene_pos_interval_tree[beg_b0 : end_b0 + 0.5] = line

					if beg_b0 not in self.gene_pos_beg_coord_dict:
						self.gene_pos_beg_coord_dict[beg_b0] = [(end_b0, line)]
					else:
						self.gene_pos_beg_coord_dict[beg_b0].append((end_b0, line))

					if end_b0 not in self.gene_pos_end_coord_dict:
						self.gene_pos_end_coord_dict[end_b0] = [(beg_b0, line)]
					else:
						self.gene_pos_end_coord_dict[end_b0].append((beg_b0, line))

				if strand == '-':
					self.gene_neg_interval_tree[beg_b0 : end_b0 + 0.5] = line

					neg_beg_b0 = end_b0
					neg_end_b0 = beg_b0

					if neg_beg_b0 not in self.gene_neg_beg_coord_dict:
						self.gene_neg_beg_coord_dict[neg_beg_b0] = [(neg_end_b0, line)]
					else:
						self.gene_neg_beg_coord_dict[neg_beg_b0].append((neg_end_b0, line))

					if neg_end_b0 not in self.gene_neg_end_coord_dict:
						self.gene_neg_end_coord_dict[neg_end_b0] = [(neg_beg_b0, line)]
					else:
						self.gene_neg_end_coord_dict[neg_end_b0].append((neg_beg_b0, line))
		in_gff_file.close()

		pos_gene_beg_sorted = sorted(self.gene_pos_beg_coord_dict.keys())
		neg_gene_beg_sorted = sorted(self.gene_neg_beg_coord_dict.keys())

		self.first_pos_gene = (pos_gene_beg_sorted[0], 
							   self.gene_pos_beg_coord_dict[pos_gene_beg_sorted[0]])
		self.last_pos_gene = (pos_gene_beg_sorted[-1], 
							  self.gene_pos_beg_coord_dict[pos_gene_beg_sorted[-1]])

		self.first_neg_gene = (neg_gene_beg_sorted[0], 
							   self.gene_neg_beg_coord_dict[neg_gene_beg_sorted[0]])
		self.last_neg_gene = (neg_gene_beg_sorted[-1], 
							  self.gene_neg_beg_coord_dict[neg_gene_beg_sorted[-1]])

	def get_overlap_trna_list(self, coord_b0, strand):
		overlap_trna_set = self.trna_interval_tree[coord_b0]
		overlap_trna_list = []
		for trna_interval in overlap_trna_set:
			if strand == trna_interval[2][6]:
				overlap_trna_list.append(trna_interval)
		return overlap_trna_list

	def get_nearest_upstream_gene(self, coord_b0, strand):
		if strand == '+':
			for scan_coord_b0 in xrange(coord_b0, -1, -1):
				if scan_coord_b0 in self.gene_pos_beg_coord_dict:
					gene_end_tuple_list = self.gene_pos_beg_coord_dict[scan_coord_b0]
					filtered_gene_end_tuple_list = gene_end_tuple_list

					distance = scan_coord_b0 - coord_b0
					return (distance, filtered_gene_end_tuple_list)

			distance = -(self.ref_length - self.last_pos_gene[0] + coord_b0)
			return (distance, self.last_pos_gene[1])
		elif strand == '-':
			for scan_coord_b0 in xrange(coord_b0, self.ref_length):
				if scan_coord_b0 in self.gene_neg_beg_coord_dict:
					neg_gene_end_tuple_list = self.gene_neg_beg_coord_dict[scan_coord_b0]
					filtered_neg_gene_end_tuple_list = neg_gene_end_tuple_list

					distance = scan_coord_b0 - coord_b0
					return (distance, filtered_neg_gene_end_tuple_list)
			distance = self.first_neg_gene[0] + self.ref_length - coord_b0
			return (distance, self.first_neg_gene[1])

		raise ValueError("Unexpected Strand %s" % strand)

	def get_nearest_downstream_gene(self, coord_b0, strand):
		if strand == '+':
			for scan_coord_b0 in xrange(coord_b0, self.ref_length):
				if scan_coord_b0 in self.gene_pos_beg_coord_dict:
					gene_end_tuple_list = self.gene_pos_beg_coord_dict[scan_coord_b0]
					distance = scan_coord_b0 - coord_b0
					return (distance, gene_end_tuple_list)

			distance = self.first_pos_gene[0] + self.ref_length - coord_b0
			return(distance, self.first_pos_gene[1])
		elif strand == '-':
			for scan_coord_b0 in xrange(coord_b0, -1, -1):
				if scan_coord_b0 in self.gene_neg_beg_coord_dict:
					gene_end_tuple_list = self.gene_neg_beg_coord_dict[scan_coord_b0]
					distance = scan_coord_b0 - coord_b0
					return (distance, gene_end_tuple_list)
			distance = -(self.ref_length - self.last_neg_gene[0] + coord_b0)
			return (distance, self.last_neg_gene[1])

		raise ValueError("Unexpected Strand %s" % strand)






		
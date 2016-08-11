import collections

# aln_lcoord_b0: 0 based coordinate of the alignment start (left most bp of 
# alignment)
# num_alt_aln: number of alternative alignment positions of the same read
# mm_desc: mismatch description (comma sep)
BowtieRecord = collections.namedtuple('BowtieRecord', 
                                      ['seqid', 'strand', 'refid', 
                                       'aln_lcoord_b0', 'seq', 'qscore', 
                                       'num_alt_aln', 'mm_desc'])

def iterate_bowtie_out_file(bt_fn):
    bt_file = open(bt_fn, 'r')

    for line in bt_file:
        fields = line.strip('\n').split('\t')
        if len(fields) != 8:
            raise ValueError("Number of fields not equal to 8: %s" % line)
        rec = BowtieRecord(fields[0], fields[1], fields[2], int(fields[3]), 
                           fields[4], fields[5], int(fields[6]), fields[7])
        if rec.strand not in ('+', '-'):
            raise ValueError("Strand not +/-: %s" % line)
        yield rec

    bt_file.close()

# Treat genome as circular
class BowtieRecordCounter(object):
    """docstring for BowtieRecordCounter"""
    def __init__(self, ref_length):
        super(BowtieRecordCounter, self).__init__()
        self.ref_length = ref_length
        self.align_count_dict = {}

    def insert_bt_rec_ntup(self, bt_rec_ntup):
        if bt_rec_ntup.aln_lcoord_b0 >= self.ref_length:
            raise ValueError("Alignment start >= ref length. %s" % bt_rec_ntup._asdict())
        
        if bt_rec_ntup.strand == '+':
            tx_start_pos_b0 = bt_rec_ntup.aln_lcoord_b0
        else:
            tx_start_pos_b0 = bt_rec_ntup.aln_lcoord_b0 + len(bt_rec_ntup.seq) - 1
            if tx_start_pos_b0 >= self.ref_length:
                tx_start_pos_b0 -= self.ref_length

        if tx_start_pos_b0 not in self.align_count_dict:
            self.align_count_dict[tx_start_pos_b0] = {'+' : 0, '-' : 0}

        self.align_count_dict[tx_start_pos_b0][bt_rec_ntup.strand] += 1

    def output_count_table(self, output_fn):
        output_file = open(output_fn, 'w')

        for key in sorted(self.align_count_dict.keys()):
            output_file.write("%d\t%d\t%d\n" % (key, 
                                                self.align_count_dict[key]['+'], 
                                                self.align_count_dict[key]['-']))

        output_file.close()
        
        
        
        
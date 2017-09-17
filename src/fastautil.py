# python2.7 notice: use of integer division. 
import dnautil
import collections as cl

FastaQuery = cl.namedtuple('FastaQuery', ['coord', 'num_up_bp', 'num_down_bp', 'strand'])
FastaQueryResult = cl.namedtuple('FastaQuery', ['coord', 'num_up_bp', 'num_down_bp', 'strand', 'seq'])

class SingleEntryFastaFile(object):
    """Simple class for manipulating single record fasta file"""
    def __init__(self, file_path):
        super(SingleEntryFastaFile, self).__init__()
        self.file_path = file_path
        header, ref_seq = self.get_fasta_header_seq(file_path)
        self.header = header
        self.ref_seq = ref_seq
        self.ref_seq_len = len(ref_seq)

    def get_seq(self, start_ind, end_ind):
        return self.ref_seq[start_ind:end_ind]

    def get_circularized_seq(self, start_ind, num_bp):
        if start_ind < 0 or start_ind >= self.ref_seq_len:
            raise ValueError("start_ind should >= 0: %s" % start_ind)
            
        if num_bp < 0:
            raise ValueError("num_bp should >= 0: %s" % num_bp)

        end_ind = start_ind + num_bp
        if end_ind <= len(self.ref_seq):
            # queried sequence
            q_seq = self.ref_seq[start_ind:end_ind]
        else:
            # when num_bp covers more than 2 copies of whole sequence
            extra_bp_after_one_complete_seq = end_ind - len(self.ref_seq)
            num_whole_seq = extra_bp_after_one_complete_seq / len(self.ref_seq)
            num_trail_bp = extra_bp_after_one_complete_seq % len(self.ref_seq)
            q_seq = self.ref_seq[start_ind:] + self.ref_seq * num_whole_seq + self.ref_seq[:num_trail_bp]

        assert len(q_seq) == num_bp
        return q_seq
        
    def get_stranded_circularized_seq(self, coord, up_bp, down_bp, strand):
        assert coord >= 0 and coord < self.ref_seq_len
        assert up_bp >= 0
        assert down_bp >= 0
        
        if strand == '+':
            left_bp = up_bp
            right_bp = down_bp
        elif strand == '-':
            left_bp = down_bp
            right_bp = up_bp
        else:
            raise ValueError("Unknown strand type: %s" % strand)

        # coord nt is the + strand nt
        # + strand right sequence can directly be extracted using self.get_circularized_seq()
        c_right_seq = self.get_circularized_seq(coord, right_bp + 1)
        c_nt = c_right_seq[0]
        right_seq = c_right_seq[1:]

        # + strand left seq
        # left seq start index might be < 0
        # abs_left_ind is the number of extra bp needed before the 
        # reference start
        # ACGTACG
        # 0123456
        # coord = 1
        # left_bp = 10
        # left_ind = -9, need 9 bp before the reference start
        # num_end_to_start_bp = abs(left_ind) % ref_seq_len = 9 % 7 = 2
        # This is the number of bp needed from the ref end to start.
        #
        # 0 <= num_end_to_start_bp <= ref_seq_len - 1
        # ACGTACT
        # 0654321
        left_ind = coord - left_bp
        if left_ind < 0:
            abs_left_ind = abs(left_ind)
            num_end_to_start_bp = abs_left_ind % self.ref_seq_len
            if num_end_to_start_bp != 0:
                left_ind = self.ref_seq_len - num_end_to_start_bp
            else:
                left_ind = 0
        
        left_seq = self.get_circularized_seq(left_ind, left_bp)

        lcr_seq = left_seq.upper() + c_nt.lower() + right_seq

        if strand == '+':
            return lcr_seq
        elif strand == '-':
            return  dnautil.reverse_complement(lcr_seq)
        else:
            raise ValueError("Unknown strand type: %s" % strand)


    @staticmethod
    def get_fasta_header_seq(file_path):
        ref_genome_file = open(file_path)
        linenum = 0
        description = ''
        ref_seq = ''
        
        for line in ref_genome_file:
            linenum += 1
            if linenum == 1:
                description = line.strip()
            else:
                if line[0] == '>':
                    raise ValueError("Input file has more than 1 record (>xxx).")
                ref_seq += line.strip()

        ref_genome_file.close()

        return (description, ref_seq)

        
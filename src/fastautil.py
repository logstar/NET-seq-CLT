# python2.7 notice: use of integer division. 
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
        if num_bp < 0:
            raise ValueError("num_bp must be positive")
        end_ind = start_ind + num_bp
        if end_ind <= len(self.ref_seq):
            return self.ref_seq[start_ind:end_ind]
        else:
            extra_bp_after_one_complete_seq = end_ind - len(self.ref_seq)
            num_whole_seq = extra_bp_after_one_complete_seq / len(self.ref_seq)
            num_trail_bp = extra_bp_after_one_complete_seq % len(self.ref_seq)
            return self.ref_seq[start_ind:] + self.ref_seq * num_whole_seq + self.ref_seq[:num_trail_bp]

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

        
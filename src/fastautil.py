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

    @staticmethod
    def get_fasta_header_seq(file_path):
        ref_genome_file = open(file_path)
        linenum = 0
        description = ''
        ref_seq = ''
        
        for line in ref_genome_file:
            linenum += 1
            if linenum == 1:
                description = line[:-1]
            else:
                if line[0] == '>':
                    raise ValueError("Input file has more than 1 record (>xxx).")
                ref_seq += line[:-1]

        ref_genome_file.close()

        ref_len = len(ref_seq)

        return (description, ref_seq)

        
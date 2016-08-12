import collections
import fastautil

ReadCountTuple = collections.namedtuple('ReadCountTuple', 
                                        ['fwd_cnt', 'rev_cnt'])

# Forward reference genome coordinate count table.
# Tab delimited three fields
# Coordinate: 0-based
# number of reads started at + strand
# number of reads started at - strand
class CoordCountTable(object):
    """docstring for CoordinateCountTable"""
    def __init__(self, tbl_fn, ref_fn):
        super(CoordCountTable, self).__init__()
        self.tbl_fn = tbl_fn
        self.ref_fn = ref_fn

        self.ref_genome = fastautil.SingleEntryFastaFile(ref_fn)

        self.coord_cnt_dict = {}
        with open(tbl_fn, 'r') as tbl_file:
            for line in tbl_file:
                fields = line.strip('\n').split('\t')
                coord = int(fields[0])
                fwd_cnt = int(fields[1])
                rev_cnt = int(fields[2])

                assert coord not in self.coord_cnt_dict, 'Duplicate coordinates. ' + line
                assert coord >= 0, 'Coordinate < 0. ' + line
                assert coord < self.ref_genome.ref_seq_len, 'Coordinate >= ref genome len. ' + line

                assert fwd_cnt >= 0, 'Forward count < 0. ' + line
                assert rev_cnt >= 0, 'Reverse count < 0. ' + line

                self.coord_cnt_dict[coord] = ReadCountTuple(fwd_cnt = fwd_cnt, 
                                                            rev_cnt = rev_cnt)

    # Coordinate is 0-based
    def get_up_down_stream_seq_cnt_tuple(self, coord, strand, num_up_bp, num_down_bp):
        left_range, right_range = self.get_left_right_range(coord, strand, 
                                                            num_up_bp, 
                                                            num_down_bp)

        report_range = left_range + [coord] + right_range
        
        if strand == '+':
            count_list = [self.coord_cnt_dict[x].fwd_cnt if x in self.coord_cnt_dict else 0 
                          for x in report_range]
            seq = ''.join([self.ref_genome.ref_seq[x].upper() if x != coord 
                           else self.ref_genome.ref_seq[x].lower() for x in report_range])
        else:
            count_list = [self.coord_cnt_dict[x].rev_cnt if x in self.coord_cnt_dict else 0 
                          for x in report_range][::-1]
            seq = self.rev_comp(''.join([self.ref_genome.ref_seq[x].upper() if x != coord 
                           else self.ref_genome.ref_seq[x].lower() for x in report_range]))

        return (seq, tuple(count_list))
        
    def get_left_right_range(self, coord, strand, num_up_bp, num_down_bp):
        assert coord < self.ref_genome.ref_seq_len, 'Coord %s >= ref genome len' % coord
        assert coord >= 0, 'Coord %s < 0' % coord
        
        assert strand in ('+', '-'), 'Strand %s not equal to + or -' % strand
        
        assert num_up_bp >= 0, '# upstream bp %s < 0' % num_up_bp
        assert num_down_bp >= 0, '# downstream bp %s < 0' % num_down_bp

        # inclusive coord on forward strand
        if strand == '+':
            left_coord = coord - num_up_bp
            right_coord = coord + num_down_bp
        else:
            left_coord = coord - num_down_bp
            right_coord = coord + num_up_bp

        if left_coord < 0:
            left_range = range(self.ref_genome.ref_seq_len + left_coord, 
                self.ref_genome.ref_seq_len, 1) + range(0, coord, 1)
        else:
            left_range = range(left_coord, coord, 1)

        if right_coord > self.ref_genome.ref_seq_len:
            right_range = (range(coord + 1, self.ref_genome.ref_seq_len, 1) + 
                range(0, right_coord - self.ref_genome.ref_seq_len + 1, 1))
        else:
            right_range = range(coord + 1, right_coord + 1, 1)

        if strand == '+':
            assert len(left_range) == num_up_bp
            assert len(right_range) == num_down_bp
        else:
            assert len(left_range) == num_down_bp
            assert len(right_range) == num_up_bp
        
        return (left_range, right_range)
        
    @staticmethod
    def rev_comp(seq):
        refd = dict(zip('ACGTacgt','TGCAtgca'))
        cseq = ''.join(map(lambda i: refd[i], seq))
        rcseq = cseq[::-1]
        return rcseq        



        
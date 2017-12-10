# python2.7 only
import collections

# check match from left to right.
def get_seq_match_len(seqx, seqy):
    match_len = 0
    for i in xrange(min(len(seqx), len(seqy))):
        if seqx[i] == seqy[i]:
            match_len += 1
        else:
            break
    
    return match_len

class TargetSeqMatchLenCounter(object):
    """docstring for SeqMatchLenCounter"""
    def __init__(self, target_seq):
        super(TargetSeqMatchLenCounter, self).__init__()
        self._target_seq = target_seq
        self._mlen_cnt_dict = collections.defaultdict(lambda: 0)

    def add_query_seq(self, query_seq):
        mlen = get_seq_match_len(self._target_seq, query_seq)
        self._mlen_cnt_dict[mlen] += 1

    # max_mlen: max output matched length, default 50.
    # min_mlen: min output matched length, default 10.
    def to_list(self, min_mlen = 10, max_mlen = 50):
        if min_mlen < 0 or max_mlen < 0:
            raise ValueError("min_mlen and max_mlen should both >= 0")

        if min_mlen >= max_mlen:
            raise ValueError("min_mlen should < max_mlen")

        tseq_cnt_list = [self._target_seq] + map(lambda mlen: (mlen, self._mlen_cnt_dict[mlen]),
                                                 xrange(min_mlen, max_mlen + 1))

        return tseq_cnt_list


# Check from beginning to the end of the target seq,
# whether the querying sequence is matched or not.
# Count the number of matches. 
class TargetSeqMatchCounter(object):
    """docstring for SeqMatchLenCounter"""
    def __init__(self, target_seq):
        super(TargetSeqMatchCounter, self).__init__()
        self._target_seq = target_seq
        self._target_seq_len = len(target_seq)
        # { match_start_ind : match_count }
        self._m_sind_cnt_dict = collections.defaultdict(lambda: 0)

    def add_query_seq(self, query_seq):
        if query_seq == '':
            return
        
        for start_ind in xrange(0, self._target_seq_len - len(query_seq) + 1):
            if self._target_seq[start_ind : start_ind + len(query_seq)] == query_seq:
                self._m_sind_cnt_dict[start_ind] += 1

    # Because the match checking requires >= n bases match, so certain length of the 
    # sequence end should all have 0 counts. Skip outputting them to prevent confusion. 
    def to_list(self, end_skip_len = 0):
        if end_skip_len < 0:
            raise ValueError("end_skip_len should >= 0")

        if end_skip_len >= self._target_seq_len:
            raise ValueError("end_skip_len should < target seq len.")

        tseq_cnt_list = [self._target_seq] + map(lambda sind: (sind, self._m_sind_cnt_dict[sind]),
                                                 xrange(0, self._target_seq_len - end_skip_len))
        return tseq_cnt_list


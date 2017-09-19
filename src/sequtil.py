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

    def to_list(self):
        tseq_cnt_list = [self._target_seq] + map(lambda i: self._mlen_cnt_dict[i],
                                                 xrange(len(self._target_seq) + 1))

        return tseq_cnt_list
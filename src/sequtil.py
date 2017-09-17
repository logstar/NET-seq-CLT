# python2.7 only

# check match from left to right.
def get_seq_match_len(seqx, seqy):
    match_len = 0
    for i in xrange(min(len(seqx), len(seqy))):
        if seqx[i] == seqy[i]:
            match_len += 1
        else:
            break
    
    return match_len

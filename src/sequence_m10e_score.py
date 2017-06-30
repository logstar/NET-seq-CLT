#!/usr/bin/env python2.7
import sys
import fastautil
import dnautil

# Assume that the genome is circular.

# Get the m10e score for all possible contiguous sequence. 
# In order to change the m10e sequence, following code needs to be checked:
# 1. indices for checking in m10e_score()
# 2. window range in window_m10e_score()

# Calculate the m10e score of a 6bp sequence.
def m10e_score(seq):
    m10e = 'TATAAT'
    if len(seq) != len(m10e):
        raise ValueError('Sequence should have the same length as m10e: %s' % seq)

    # hardcode here for initial testing and speed
    if seq[1] == m10e[1] and seq[5] == m10e[5]:
        seq_m10e_score = 2
        for i in (0, 2, 3, 4):
            if seq[i] == m10e[i]:
                seq_m10e_score += 1
        return seq_m10e_score
    else:
        return 0

# Always treat sequence as linear
# Calculate the m10e score of a sequence from the beginning to the end.
# Treat the genome as circular.
# Return: 6bp window score vector. Firt score is the window of the first bp. 
def consecutive_m10e_score(seq):
    m10e_len = 6
    consecutive_m10e_score_list = []
    for i in xrange(len(seq) - m10e_len + 1):
        consecutive_m10e_score_list.append(m10e_score(seq[i:i+m10e_len]))
    return consecutive_m10e_score_list

# Always treat sequence as linear
# 9bp. 4 scores. Frame 7bp. 3 frames. 
# XXXXXXXXX
#    ------
# [     ]
#  [     ]
#   [     ]
# Output a table of all frame m10e window scores
# pos
# frame seq
# scores
def frame_fmt_m10e_score(seq, frame_size):
    m10e_len = 6
    if frame_size < m10e_len:
        raise ValueError('Frame size should >= 6. ')
    consecutive_m10e_score_list = consecutive_m10e_score(seq)
    frame_seq_m10e_score_list = []
    for i in xrange(len(seq) - frame_size + 1):
        frame_seq = seq[i : i + frame_size]
        frame_score = tuple(consecutive_m10e_score_list[i : i + frame_size - m10e_len + 1])
        frame_seq_m10e_score_list.append((i, frame_seq, frame_score))
    return frame_seq_m10e_score_list


def is_valid_seq_type(seq_type):
    if seq_type in 'cl':
        return True
    else:
        return False


def format_frame_m10e_score_tup(tup, strand):
    fmted_string = str(tup[0]) + '\t' + strand + '\t' + tup[1] + '\t' + '\t'.join(map(str, tup[2]))
    return fmted_string


# In order to handle circular seuqnce, simply append n bp beginning sequence at the end. 
# n is the frame size - 1, so that the last bp will also have a frame.
# Tab-delimited output format:
# 1. pos 0-based
# 2. strand
# 3. frame sequence
# 4-. score of the window
def main():
    argv = sys.argv
    hlp_msg = 'Usage:\n%s <frame size> <circular/linear> <single entry fasta file> \
<output position score file>\n\
<carcular/linear>: c for circular. l for linear.\n' % argv[0]

    if len(argv) != 5:
        sys.stderr.write(hlp_msg)
        return -1

    frame_size = int(argv[1])
    seq_type = argv[2]

    if not is_valid_seq_type(seq_type):
        sys.stderr.write(hlp_msg)
        return -2

    ifn = argv[3]
    ifa_file = fastautil.SingleEntryFastaFile(ifn)

    if seq_type == 'c':
        ifa_seq = ifa_file.get_circularized_seq(0, len(ifa_file.ref_seq) + frame_size - 1)
    elif seq_type == 'l':
        ifa_seq = ifa_file.ref_seq
    else:
        raise ValueError("Invalid sequence type: %s" % seq_type)

    ifa_seq_comp = dnautil.complement(ifa_seq)
    ifa_seq_frame_m10e_score_list = frame_fmt_m10e_score(ifa_seq, 36)
    ifa_seq_comp_frame_m10e_score_list = frame_fmt_m10e_score(ifa_seq_comp, 36)
    
    ofn = argv[4]
    with open(ofn, 'w') as ofile:
        # tup: (pos:Int, seq:String, score:IntTuple)
        for frm_scr_tup in ifa_seq_frame_m10e_score_list:
            ofile.write(format_frame_m10e_score_tup(frm_scr_tup, '+') + '\n')

        for frm_scr_tup in ifa_seq_comp_frame_m10e_score_list:
            ofile.write(format_frame_m10e_score_tup(frm_scr_tup, '-') + '\n')


if __name__ == "__main__":
    main()

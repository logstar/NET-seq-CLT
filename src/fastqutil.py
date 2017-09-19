# Input file name of fastq file, yield every 4 lines
# Yield a 4-element tuple (seqid, seq, rseqid, qscore)
# Lines are stripped
def iterate_fastq_file(seq_fn):
    with open(seq_fn, 'r') as seq_file:
        parsed_line_num = 0
        seqid = ''
        seq = ''
        rseqid = ''
        qscore = ''
        for line in seq_file:
            parsed_line_num += 1

            ln_mod_4 = parsed_line_num % 4

            if ln_mod_4 == 1:
                seqid = line.strip()
            elif ln_mod_4 == 2:
                seq = line.strip()
            elif ln_mod_4 == 3:
                rseqid = line.strip()
            elif ln_mod_4 == 0:
                qscore = line.strip()
                yield (seqid, seq, rseqid, qscore)
            else:
                seq_file.close()
                raise ValueError("Should be impossible: parsed_line_num % 4 not in \
    [1, 2, 3, 0]\n")
                break

    if parsed_line_num % 4 != 0:
        raise ValueError("Number of line is not a multiple of 4")


def int_to_qscore(qint):
    if not isinstance(qint, int):
        raise ValueError('Input must be an integer.')
        
    if qint < 0 or qint > 93:
        raise ValueError('Quality score cutoff should >= 0 and <= 93')

    qchar = chr(qint + 33)
    return(qchar)

#!/usr/bin/env python2.7
import argparse
import fastautil


def get_fasta_query_list(q_fn):
    fasta_query_list = []
    with open(q_fn, 'r') as q_file:
        for line in q_file:
            fields = line.strip().split()
            assert len(fields) == 4
            coord = int(fields[0])
            num_up_bp = int(fields[1])
            num_down_bp = int(fields[2])
            strand = fields[3]
            assert strand in '+-'
            fasta_query_list.append(fastautil.FastaQuery(coord = coord, num_up_bp = num_up_bp,
                                                         num_down_bp = num_down_bp, strand = strand))

    return fasta_query_list

def query_fasta(se_fa_file, fasta_query):
    queried_seq = se_fa_file.get_stranded_circularized_seq(coord = fasta_query.coord,
        up_bp = fasta_query.num_up_bp, down_bp = fasta_query.num_down_bp, strand = fasta_query.strand)

    query_result = fastautil.FastaQueryResult(coord = fasta_query.coord, 
                                              num_up_bp = fasta_query.num_up_bp, 
                                              num_down_bp = fasta_query.num_down_bp, 
                                              strand = fasta_query.strand, 
                                              seq = queried_seq)
    return query_result

def main():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('se_fa_file', metavar = '<single entry fasta file>')
    arg_parser.add_argument('q_file_name', metavar = '<query coord numbp strand file>',
                            help = 'A tab delimited file listing all queries. ' 
                                   'Fields from left to right are: '
                                   'coordinate, number of upstream bp, '
                                   'number of downstream bp, strand. No header.')
    arg_parser.add_argument('o_file_name', metavar = '<output file>')

    args = arg_parser.parse_args()

    se_fa_file = fastautil.SingleEntryFastaFile(args.se_fa_file)
    fasta_query_list = get_fasta_query_list(args.q_file_name)

    with open(args.o_file_name, 'w') as ofile:
        for fasta_query in fasta_query_list:
            ofile.write('\t'.join(map(str, query_fasta(se_fa_file, fasta_query))) + '\n')

    return

if __name__ == '__main__':
    main()
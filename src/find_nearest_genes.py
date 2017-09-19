#!/usr/bin/env python2.7
import sys
from gffutils import GFF
import argparse

def find_nearest_genes(in_tss_fn, o_tss_fn, ref_gff):
    in_tss_file = open(in_tss_fn, 'r')
    o_tss_file = open(o_tss_fn, 'w')

    for line in in_tss_file:
        fields = line[:-1].split('\t')
        coord = int(fields[0])
        strand = fields[1]
        nearest_upstream_gene_tup = ref_gff.get_nearest_upstream_gene(coord, strand)
        nearest_downstream_gene_tup = ref_gff.get_nearest_downstream_gene(coord, strand)

        o_tss_file.write(line.strip())

        o_tss_file.write("\t%d\t" % nearest_upstream_gene_tup[0])
        o_tss_file.write('|'.join(map(lambda tup: tup[1][:-1].replace('\t', ' '), nearest_upstream_gene_tup[1])))
        o_tss_file.write("\t%d\t" % nearest_downstream_gene_tup[0])
        o_tss_file.write('|'.join(map(lambda tup: tup[1][:-1].replace('\t', ' '), nearest_downstream_gene_tup[1])))
        o_tss_file.write('\n')
    in_tss_file.close()
    o_tss_file.close()

def main():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('in_tss_fn', metavar = '<input TSS file (0-based)>',
                             help = 'TSS file has at least two fields, which are '
                                    '0-based coordinate and strand from left '
                                    'to right.')

    arg_parser.add_argument('in_gff_fn', metavar = '<input reference GFF file>',
    	                    help = 'Gene annotation file in GFF3 format. ')

    arg_parser.add_argument('o_tss_ng_fn', metavar = '<output TSS nearest gene file>',
    	                    help = 'Output file including the nearest upstream and '
    	                           'downstream genes appended at the end of each line.')

    args = arg_parser.parse_args()

    ref_gff = GFF(args.in_gff_fn)

    find_nearest_genes(args.in_tss_fn, args.o_tss_ng_fn, ref_gff)



if __name__ == "__main__":
	main()
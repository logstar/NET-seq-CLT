#!/usr/bin/env python2.7
import sys
from gffutils import GFF
def find_opposite_strand_nearest_genes(in_tss_fn, o_tss_fn, ref_gff):
    in_tss_file = open(in_tss_fn, 'r')
    o_tss_file = open(o_tss_fn, 'w')

    for line in in_tss_file:
        fields = line[:-1].split('\t')
        coord = int(fields[0])
        strand = fields[1]
        if strand == '+':
            opposite_strand = '-'
        elif strand == '-':
            opposite_strand = '+'
        else:
            raise ValueError("Unknown strand for %s" % line)
        
        nearest_upstream_gene_tup = ref_gff.get_nearest_upstream_gene(coord, opposite_strand)
        nearest_downstream_gene_tup = ref_gff.get_nearest_downstream_gene(coord, opposite_strand)

        o_tss_file.write(line.strip())

        o_tss_file.write("\t%d\t" % nearest_upstream_gene_tup[0])
        o_tss_file.write('|'.join(map(lambda tup: tup[1][:-1].replace('\t', ' '), nearest_upstream_gene_tup[1])))
        o_tss_file.write("\t%d\t" % nearest_downstream_gene_tup[0])
        o_tss_file.write('|'.join(map(lambda tup: tup[1][:-1].replace('\t', ' '), nearest_downstream_gene_tup[1])))
        o_tss_file.write('\n')
    in_tss_file.close()
    o_tss_file.close()

if __name__ == "__main__":
    argv = sys.argv
    if len(argv) != 4:
        sys.stderr.write("Usage:\n%s <input TSS file (0 based)>\
 <input reference GFF file> <output TSS nearest gene file> \n" % argv[0])
        sys.exit(2)

    in_tss_fn = argv[1]
    in_gff_fn = argv[2]
    o_tss_fn = argv[3]

    ref_gff = GFF(in_gff_fn)

    find_opposite_strand_nearest_genes(in_tss_fn, o_tss_fn, ref_gff)


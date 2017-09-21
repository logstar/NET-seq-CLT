#!/usr/bin/env python2.7
import argparse

def gff_remove_pseudo_gene(in_gff_fn, out_gff_fn):
    with open(in_gff_fn, 'r') as ifile:
        with open(out_gff_fn, 'w') as ofile:
            for line in ifile:
                if 'pseudo=true' not in line:
                    ofile.write(line)

    return

def main():
    arg_parser = argparse.ArgumentParser(description = 'Simply remove the lines with \'pseudo=true\'.')

    arg_parser.add_argument('in_gff_fn', metavar = '<input GFF3 file with pseudo genes>')
    arg_parser.add_argument('out_gff_fn', metavar = '<out GFF3 file without pseudo genes>')

    args = arg_parser.parse_args()
    
    gff_remove_pseudo_gene(args.in_gff_fn, args.out_gff_fn)

if __name__ == '__main__':
    main()
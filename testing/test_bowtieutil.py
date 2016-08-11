#!/usr/bin/env python2.7
import sys
import unittest
import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

sys.path.append('../src')
import bowtieutil as bu

class TestIterator(unittest.TestCase):
    def test_tuple_order(self):
        t_btfn = 'data/test_bu1.txt'
        with open(t_btfn, 'w') as t_btfile:
            t_btfile.write(self.get_bowtie_out_record('1','+','3','4','5','6','7','8'))
            t_btfile.write(self.get_bowtie_out_record('2','+','3','4','5','6','7','8'))
            t_btfile.write(self.get_bowtie_out_record('3','+','3','4','5','6','7','8'))

        linenum = 0
        for rec in bu.iterate_bowtie_out_file(t_btfn):
            linenum += 1
            self.assertEqual(dict(rec._asdict()), 
                             dict(zip(('seqid', 'strand', 'refid', 
                                       'aln_lcoord_b0', 'seq', 'qscore', 
                                       'num_alt_aln', 'mm_desc'), 
                                      (str(linenum), '+','3','4','5','6','7','8'))))

        with open(t_btfn, 'w') as t_btfile:
            t_btfile.write(self.get_bowtie_out_record('1','+','3','4','5','6','',''))
            t_btfile.write(self.get_bowtie_out_record('2','+','3','4','5','6','',''))
            t_btfile.write(self.get_bowtie_out_record('3','+','3','4','5','6','',''))

        linenum = 0
        for rec in bu.iterate_bowtie_out_file(t_btfn):
            linenum += 1
            self.assertEqual(dict(rec._asdict()), 
                             dict(zip(('seqid', 'strand', 'refid', 
                                       'aln_lcoord_b0', 'seq', 'qscore', 
                                       'num_alt_aln', 'mm_desc'), 
                                      (str(linenum), '+','3','4','5','6','',''))))

    def test_wrong_num_field(self):
        t_btfn = 'data/test_bu2.txt'
        with open(t_btfn, 'w') as t_btfile:
            t_btfile.write('\t'.join(('1','+','3','4','5','6','7')) + '\n')

        with self.assertRaises(ValueError):
            for rec in bu.iterate_bowtie_out_file(t_btfn):
                pass

    def test_wrong_strand(self):
        t_btfn = 'data/test_bu2.txt'
        with open(t_btfn, 'w') as t_btfile:
            t_btfile.write('\t'.join(('1','2','3','4','5','6','7', '8')) + '\n')

        with self.assertRaises(ValueError):
            for rec in bu.iterate_bowtie_out_file(t_btfn):
                pass

        

    def test_line_num(self):
        t_btfn = 'data/test_bu1.txt'
        with open(t_btfn, 'w') as t_btfile:
            t_btfile.write(self.get_bowtie_out_record('1','+','3','4','5','6','7',''))
            t_btfile.write(self.get_bowtie_out_record('2','+','3','4','5','6','7','8'))
            t_btfile.write(self.get_bowtie_out_record('3','+','3','4','5','6','7','8'))

        linenum = 0
        for rec in bu.iterate_bowtie_out_file(t_btfn):
            linenum += 1

        self.assertEqual(linenum, 3)

        with open(t_btfn, 'w') as t_btfile:
            t_btfile.write(self.get_bowtie_out_record('1','+','3','4','5','6','7','8'))

        linenum = 0
        for rec in bu.iterate_bowtie_out_file(t_btfn):
            linenum += 1

        self.assertEqual(linenum, 1)

        with open(t_btfn, 'w') as t_btfile:
            pass
        linenum = 0
        for rec in bu.iterate_bowtie_out_file(t_btfn):
            linenum += 1

        self.assertEqual(linenum, 0)

    @staticmethod
    def get_bowtie_out_record(seqid, strand, refid, aln_lcoord_b0, seq, qscore, 
                              num_alt_aln, mm_desc):
        rec = '\t'.join((seqid, strand, refid, aln_lcoord_b0, seq, qscore, 
                         num_alt_aln, mm_desc)) + '\n'
        return rec

if __name__ == '__main__':
    unittest.main()
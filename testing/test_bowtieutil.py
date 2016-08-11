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
                                      (str(linenum), '+','3',4,'5','6',7,'8'))))

        with open(t_btfn, 'w') as t_btfile:
            t_btfile.write(self.get_bowtie_out_record('1','+','3','4','5','6','0',''))
            t_btfile.write(self.get_bowtie_out_record('2','+','3','4','5','6','0',''))
            t_btfile.write(self.get_bowtie_out_record('3','+','3','4','5','6','0',''))

        linenum = 0
        for rec in bu.iterate_bowtie_out_file(t_btfn):
            linenum += 1
            self.assertEqual(dict(rec._asdict()), 
                             dict(zip(('seqid', 'strand', 'refid', 
                                       'aln_lcoord_b0', 'seq', 'qscore', 
                                       'num_alt_aln', 'mm_desc'), 
                                      (str(linenum), '+','3',4,'5','6',0,''))))

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

class TestBowtieRecordCounter(unittest.TestCase):
    def test_counting(self):
        align_counter = bu.BowtieRecordCounter(100)
        self.assertEqual(align_counter.ref_length, 100)

        # Test +
        align_counter.insert_bt_rec_ntup(bu.BowtieRecord._make(
            ('1', '+', 'tref', 0, 'ACGTACGT', chr(50) * 8, 0, '')))
        self.assertEqual(align_counter.align_count_dict[0]['+'], 1)
        align_counter.insert_bt_rec_ntup(bu.BowtieRecord._make(
            ('1', '+', 'tref', 0, 'ACGTACGT', chr(50) * 8, 0, '')))
        self.assertEqual(align_counter.align_count_dict[0]['+'], 2)
        align_counter.insert_bt_rec_ntup(bu.BowtieRecord._make(
            ('1', '+', 'tref', 0, 'ACGTACGT', chr(50) * 8, 0, '')))
        self.assertEqual(align_counter.align_count_dict[0]['+'], 3)

        # Test -
        align_counter.insert_bt_rec_ntup(bu.BowtieRecord._make(
            ('1', '-', 'tref', 0, 'ACGTACGT', chr(50) * 8, 0, '')))
        self.assertEqual(align_counter.align_count_dict[7]['-'], 1)
        align_counter.insert_bt_rec_ntup(bu.BowtieRecord._make(
            ('1', '-', 'tref', 0, 'ACGTACGT', chr(50) * 8, 0, '')))
        self.assertEqual(align_counter.align_count_dict[7]['-'], 2)
        self.assertEqual(align_counter.align_count_dict[0]['+'], 3)

        # Test circular -
        align_counter.insert_bt_rec_ntup(bu.BowtieRecord._make(
            ('1', '-', 'tref', 99, 'ACGTACGT', chr(50) * 8, 0, '')))
        self.assertEqual(align_counter.align_count_dict[6]['-'], 1)
        align_counter.insert_bt_rec_ntup(bu.BowtieRecord._make(
            ('1', '-', 'tref', 99, 'ACGTACGT', chr(50) * 8, 0, '')))
        self.assertEqual(align_counter.align_count_dict[6]['-'], 2)
        align_counter.insert_bt_rec_ntup(bu.BowtieRecord._make(
            ('1', '-', 'tref', 99, 'ACGTACGT', chr(50) * 8, 0, '')))
        self.assertEqual(align_counter.align_count_dict[6]['-'], 3)

        align_counter.insert_bt_rec_ntup(bu.BowtieRecord._make(
            ('1', '-', 'tref', 93, 'ACGTACGT', chr(50) * 8, 0, '')))
        self.assertEqual(align_counter.align_count_dict[0]['-'], 1)
        align_counter.insert_bt_rec_ntup(bu.BowtieRecord._make(
            ('1', '-', 'tref', 93, 'ACGTACGT', chr(50) * 8, 0, '')))
        self.assertEqual(align_counter.align_count_dict[0]['-'], 2)
        align_counter.insert_bt_rec_ntup(bu.BowtieRecord._make(
            ('1', '-', 'tref', 93, 'ACGTACGT', chr(50) * 8, 0, '')))
        self.assertEqual(align_counter.align_count_dict[0]['-'], 3)

    def test_output(self):
        align_counter = bu.BowtieRecordCounter(100)
        self.assertEqual(align_counter.ref_length, 100)

        # Test +
        align_counter.insert_bt_rec_ntup(bu.BowtieRecord._make(
            ('1', '+', 'tref', 0, 'ACGTACGT', chr(50) * 8, 0, '')))
        self.assertEqual(align_counter.align_count_dict[0]['+'], 1)
        align_counter.insert_bt_rec_ntup(bu.BowtieRecord._make(
            ('1', '+', 'tref', 0, 'ACGTACGT', chr(50) * 8, 0, '')))
        self.assertEqual(align_counter.align_count_dict[0]['+'], 2)

        # Test -
        align_counter.insert_bt_rec_ntup(bu.BowtieRecord._make(
            ('1', '-', 'tref', 0, 'ACGTACGT', chr(50) * 8, 0, '')))
        self.assertEqual(align_counter.align_count_dict[7]['-'], 1)

        # Test circular -
        align_counter.insert_bt_rec_ntup(bu.BowtieRecord._make(
            ('1', '-', 'tref', 99, 'ACGTACGT', chr(50) * 8, 0, '')))
        self.assertEqual(align_counter.align_count_dict[6]['-'], 1)
        align_counter.insert_bt_rec_ntup(bu.BowtieRecord._make(
            ('1', '-', 'tref', 99, 'ACGTACGT', chr(50) * 8, 0, '')))
        self.assertEqual(align_counter.align_count_dict[6]['-'], 2)
        align_counter.insert_bt_rec_ntup(bu.BowtieRecord._make(
            ('1', '-', 'tref', 99, 'ACGTACGT', chr(50) * 8, 0, '')))
        self.assertEqual(align_counter.align_count_dict[6]['-'], 3)
        align_counter.insert_bt_rec_ntup(bu.BowtieRecord._make(
            ('1', '-', 'tref', 99, 'ACGTACGT', chr(50) * 8, 0, '')))
        self.assertEqual(align_counter.align_count_dict[6]['-'], 4)
        align_counter.insert_bt_rec_ntup(bu.BowtieRecord._make(
            ('1', '-', 'tref', 99, 'ACGTACGT', chr(50) * 8, 0, '')))
        self.assertEqual(align_counter.align_count_dict[6]['-'], 5)

        align_counter.insert_bt_rec_ntup(bu.BowtieRecord._make(
            ('1', '-', 'tref', 93, 'ACGTACGT', chr(50) * 8, 0, '')))
        self.assertEqual(align_counter.align_count_dict[0]['-'], 1)
        align_counter.insert_bt_rec_ntup(bu.BowtieRecord._make(
            ('1', '-', 'tref', 93, 'ACGTACGT', chr(50) * 8, 0, '')))
        self.assertEqual(align_counter.align_count_dict[0]['-'], 2)
        align_counter.insert_bt_rec_ntup(bu.BowtieRecord._make(
            ('1', '-', 'tref', 93, 'ACGTACGT', chr(50) * 8, 0, '')))
        self.assertEqual(align_counter.align_count_dict[0]['-'], 3)

        ofn = 'data/test_bu_cnt1.txt'
        align_counter.output_count_table(ofn)
        with open(ofn, 'r') as ofile:
            cnt_str = ofile.read()

        cnt_str_ref = "%d\t%d\t%d\n" % (0, 2, 3)
        cnt_str_ref += "%d\t%d\t%d\n" % (6, 0, 5)
        cnt_str_ref += "%d\t%d\t%d\n" % (7, 0, 1)

if __name__ == '__main__':
    unittest.main()

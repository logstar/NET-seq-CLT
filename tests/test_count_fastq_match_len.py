#!/usr/bin/env python2.7
import unittest
import env
import src.count_fastq_match_len as cml
import src.sequtil as su
import src.fastautil as fau

class TestAllSameMlenRange(unittest.TestCase):
    def setUp(self):
        self.tseq_cnter_list = [su.TargetSeqMatchLenCounter('ACGTACG'),
                                su.TargetSeqMatchLenCounter('AGGTACG'),
                                su.TargetSeqMatchLenCounter('ATGTACG'),
                                su.TargetSeqMatchLenCounter('ATTGTACG'),
                                su.TargetSeqMatchLenCounter('ACCGTACG')]

    def test_all_same_mlen(self):
        self.assertTrue(cml.all_same_mlen_range(map(lambda cnter: cnter.to_list(), self.tseq_cnter_list)))
        self.assertTrue(cml.all_same_mlen_range(map(lambda cnter: cnter.to_list(), self.tseq_cnter_list[:1])))
        self.assertTrue(cml.all_same_mlen_range(map(lambda cnter: cnter.to_list(0, 100), self.tseq_cnter_list)))

    def test_diff_mlen(self):
        self.assertFalse(cml.all_same_mlen_range(
            map(lambda cnter: cnter.to_list(0, 100), self.tseq_cnter_list[:3])
            + map(lambda cnter: cnter.to_list(10, 20), self.tseq_cnter_list[3:])))

        self.assertFalse(cml.all_same_mlen_range(
            map(lambda cnter: cnter.to_list(0, 100), self.tseq_cnter_list[:4])
            + map(lambda cnter: cnter.to_list(10, 20), self.tseq_cnter_list[4:])))

        self.assertFalse(cml.all_same_mlen_range(
            map(lambda cnter: cnter.to_list(0, 100), self.tseq_cnter_list[:1])
            + map(lambda cnter: cnter.to_list(10, 20), self.tseq_cnter_list[1:])))


    def test_empty_list(self):
        self.assertTrue(cml.all_same_mlen_range([]))


# fmt_fa_qres_mlen_cnt(fa_qres_list, tseq_mlen_cnt_list)
class TestFmtFaQresMlenCnt(unittest.TestCase):
    def setUp(self):
        self.tseq_cnter_list = [su.TargetSeqMatchLenCounter('ACGTACG'),
                                su.TargetSeqMatchLenCounter('AGGTACG'),
                                su.TargetSeqMatchLenCounter('ATGTACG'),
                                su.TargetSeqMatchLenCounter('ATTGTACG'),
                                su.TargetSeqMatchLenCounter('ACCGTACG')]
        
        self.tseq_same_mlen_cnt_list = map(lambda cnter: cnter.to_list(0, 5), self.tseq_cnter_list)

        self.tseq_diff_mlen_cnt_list = (map(lambda cnter: cnter.to_list(0, 100), self.tseq_cnter_list[:1]) 
            + map(lambda cnter: cnter.to_list(10, 20), self.tseq_cnter_list[1:]))

        self.fa_qres_list = [fau.FastaQueryResult(0, 3, 3, '+', 'ACGtACG'),
                             fau.FastaQueryResult(10, 3, 3, '-', 'AGGtACG'),
                             fau.FastaQueryResult(20, 3, 3, '-', 'ATGtACG'),
                             fau.FastaQueryResult(30, 4, 3, '+', 'ATTGtACG'),
                             fau.FastaQueryResult(40, 4, 3, '+', 'ACCGtACG'),
                             fau.FastaQueryResult(50, 3, 3, '-', 'ATGtACG'),
                             fau.FastaQueryResult(60, 4, 3, '+', 'ACCGtACG'),
                             fau.FastaQueryResult(70, 4, 3, '+', 'ACCGtACG')]

        # seqs not included in tseq_cnter_list
        self.fa_qres_non_exist_list = [fau.FastaQueryResult(0, 3, 3, '+', 'CCCtGGG'),
                                       fau.FastaQueryResult(0, 5, 5, '+', 'AAACCtACGAA')]

                
    def test_empty_fa_qres_list(self):
        with self.assertRaises(ValueError):
            cml.fmt_fa_qres_mlen_cnt([], self.tseq_same_mlen_cnt_list)

    def test_empty_tseq_mlen_cnt_list(self):
        with self.assertRaises(ValueError):
            cml.fmt_fa_qres_mlen_cnt(self.fa_qres_list, [])

    def test_empty_both_input(self):
        with self.assertRaises(ValueError):
            cml.fmt_fa_qres_mlen_cnt([], [])

    def test_diff_mlen(self):
        with self.assertRaises(AssertionError):
            cml.fmt_fa_qres_mlen_cnt(self.fa_qres_list, self.tseq_diff_mlen_cnt_list)

    def test_dup_tseq(self):
        with self.assertRaises(AssertionError):
            cml.fmt_fa_qres_mlen_cnt(self.fa_qres_list, 
                self.tseq_same_mlen_cnt_list + self.tseq_same_mlen_cnt_list[:2])

    def test_fa_qres_seq_not_in_tseq(self):
        with self.assertRaises(ValueError):
            cml.fmt_fa_qres_mlen_cnt(self.fa_qres_non_exist_list, self.tseq_same_mlen_cnt_list)

        with self.assertRaises(ValueError):
            cml.fmt_fa_qres_mlen_cnt(self.fa_qres_non_exist_list + self.fa_qres_list, 
                self.tseq_same_mlen_cnt_list)

    def test_invalid_r_arg(self):
        with self.assertRaises(AssertionError):
            cml.fmt_fa_qres_mlen_cnt(self.fa_qres_list, 
                self.tseq_same_mlen_cnt_list + self.tseq_same_mlen_cnt_list[:2], 
                reverse = 0)

    def test_invalid_c_arg(self):
        with self.assertRaises(AssertionError):
            cml.fmt_fa_qres_mlen_cnt(self.fa_qres_list, 
                self.tseq_same_mlen_cnt_list + self.tseq_same_mlen_cnt_list[:2], 
                complement = 0)

    def test_invalid_rc_arg(self):
        with self.assertRaises(AssertionError):
            cml.fmt_fa_qres_mlen_cnt(self.fa_qres_list, 
                self.tseq_same_mlen_cnt_list + self.tseq_same_mlen_cnt_list[:2], 
                reverse = 0, complement = 0)
        
    def test_output_format_implicit_args(self):
        o_tbl_str = cml.fmt_fa_qres_mlen_cnt(self.fa_qres_list, self.tseq_same_mlen_cnt_list)
        ref_tbl_str = '\
coord\tstrand\tseq\t0\t1\t2\t3\t4\t5\n\
0\t+\tACGtACG\t0\t0\t0\t0\t0\t0\n\
10\t-\tAGGtACG\t0\t0\t0\t0\t0\t0\n\
20\t-\tATGtACG\t0\t0\t0\t0\t0\t0\n\
30\t+\tATTGtACG\t0\t0\t0\t0\t0\t0\n\
40\t+\tACCGtACG\t0\t0\t0\t0\t0\t0\n\
50\t-\tATGtACG\t0\t0\t0\t0\t0\t0\n\
60\t+\tACCGtACG\t0\t0\t0\t0\t0\t0\n\
70\t+\tACCGtACG\t0\t0\t0\t0\t0\t0\n'
        self.assertEqual(o_tbl_str, ref_tbl_str)

    def test_output_format_explicit_args(self):
        ref_tbl_str = '\
coord\tstrand\tseq\t0\t1\t2\t3\t4\t5\n\
0\t+\tACGtACG\t0\t0\t0\t0\t0\t0\n\
10\t-\tAGGtACG\t0\t0\t0\t0\t0\t0\n\
20\t-\tATGtACG\t0\t0\t0\t0\t0\t0\n\
30\t+\tATTGtACG\t0\t0\t0\t0\t0\t0\n\
40\t+\tACCGtACG\t0\t0\t0\t0\t0\t0\n\
50\t-\tATGtACG\t0\t0\t0\t0\t0\t0\n\
60\t+\tACCGtACG\t0\t0\t0\t0\t0\t0\n\
70\t+\tACCGtACG\t0\t0\t0\t0\t0\t0\n'
        o_tbl_str_exp = cml.fmt_fa_qres_mlen_cnt(self.fa_qres_list, self.tseq_same_mlen_cnt_list, 
                                                 reverse = False, complement = False)
        self.assertEqual(o_tbl_str_exp, ref_tbl_str)

    def test_output_format_rev_seq(self):
        ref_tbl_str_rev = '\
coord\tstrand\tseq\t0\t1\t2\t3\t4\t5\n\
0\t+\tGCAtGCA\t0\t0\t0\t0\t0\t0\n\
10\t-\tGCAtGGA\t0\t0\t0\t0\t0\t0\n\
20\t-\tGCAtGTA\t0\t0\t0\t0\t0\t0\n\
30\t+\tGCAtGTTA\t0\t0\t0\t0\t0\t0\n\
40\t+\tGCAtGCCA\t0\t0\t0\t0\t0\t0\n\
50\t-\tGCAtGTA\t0\t0\t0\t0\t0\t0\n\
60\t+\tGCAtGCCA\t0\t0\t0\t0\t0\t0\n\
70\t+\tGCAtGCCA\t0\t0\t0\t0\t0\t0\n'
        o_tbl_str_rev = cml.fmt_fa_qres_mlen_cnt(self.fa_qres_list, self.tseq_same_mlen_cnt_list, 
                                                 reverse = True, complement = False)
        self.assertEqual(o_tbl_str_rev, ref_tbl_str_rev)

    def test_output_format_comp_seq(self):
        ref_tbl_str_comp = '\
coord\tstrand\tseq\t0\t1\t2\t3\t4\t5\n\
0\t+\tTGCaTGC\t0\t0\t0\t0\t0\t0\n\
10\t-\tTCCaTGC\t0\t0\t0\t0\t0\t0\n\
20\t-\tTACaTGC\t0\t0\t0\t0\t0\t0\n\
30\t+\tTAACaTGC\t0\t0\t0\t0\t0\t0\n\
40\t+\tTGGCaTGC\t0\t0\t0\t0\t0\t0\n\
50\t-\tTACaTGC\t0\t0\t0\t0\t0\t0\n\
60\t+\tTGGCaTGC\t0\t0\t0\t0\t0\t0\n\
70\t+\tTGGCaTGC\t0\t0\t0\t0\t0\t0\n'
        o_tbl_str_comp = cml.fmt_fa_qres_mlen_cnt(self.fa_qres_list, self.tseq_same_mlen_cnt_list, 
                                                  reverse = False, complement = True)
        self.assertEqual(o_tbl_str_comp, ref_tbl_str_comp)

    def test_output_format_rev_comp_seq(self):
        ref_tbl_str_rev_comp = '\
coord\tstrand\tseq\t0\t1\t2\t3\t4\t5\n\
0\t+\tCGTaCGT\t0\t0\t0\t0\t0\t0\n\
10\t-\tCGTaCCT\t0\t0\t0\t0\t0\t0\n\
20\t-\tCGTaCAT\t0\t0\t0\t0\t0\t0\n\
30\t+\tCGTaCAAT\t0\t0\t0\t0\t0\t0\n\
40\t+\tCGTaCGGT\t0\t0\t0\t0\t0\t0\n\
50\t-\tCGTaCAT\t0\t0\t0\t0\t0\t0\n\
60\t+\tCGTaCGGT\t0\t0\t0\t0\t0\t0\n\
70\t+\tCGTaCGGT\t0\t0\t0\t0\t0\t0\n'
        o_tbl_str_rev_comp = cml.fmt_fa_qres_mlen_cnt(self.fa_qres_list, self.tseq_same_mlen_cnt_list, 
                                                      reverse = True, complement = True)
        self.assertEqual(o_tbl_str_rev_comp, ref_tbl_str_rev_comp)


if __name__ == '__main__':
    unittest.main()

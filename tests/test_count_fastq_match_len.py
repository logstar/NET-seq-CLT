#!/usr/bin/env python2.7
import unittest
import env
import src.count_fastq_match_len as cml
import src.sequtil as su

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


# def fmt_tseq_rec_mlen_cnt(tseq_rec_list, tseq_mlen_cnt_list)
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

        self.tseq_rec_list = [cml.TargetSeqRecord(0, 3, 3, '+', 'ACGtACG', 'X'),
                              cml.TargetSeqRecord(10, 3, 3, '-', 'AGGtACG', 'Y'),
                              cml.TargetSeqRecord(20, 3, 3, '-', 'ATGtACG', 'ATGtACG'),
                              cml.TargetSeqRecord(30, 4, 3, '+', 'ATTGtACG', 'ATTGtACG'),
                              cml.TargetSeqRecord(40, 4, 3, '+', 'ACCGtACG', 'ACCGtACG'),
                              cml.TargetSeqRecord(50, 3, 3, '-', 'ATGtACG', 'ATGtACG'),
                              cml.TargetSeqRecord(60, 4, 3, '+', 'ACCGtACG', 'ACCGtACG'),
                              cml.TargetSeqRecord(70, 4, 3, '+', 'ACCGtACG', 'ACCGtACG')]

        # seqs not included in tseq_cnter_list
        self.tseq_rec_non_exist_list = [cml.TargetSeqRecord(0, 3, 3, '+', 'CCCtGGG', 'CCCtGGG'),
                                       cml.TargetSeqRecord(0, 5, 5, '+', 'AAACCtACGAA', 'AAACCtACGAA')]

                
    def test_empty_tseq_rec_list(self):
        with self.assertRaises(ValueError):
            cml.fmt_tseq_rec_mlen_cnt([], self.tseq_same_mlen_cnt_list)

    def test_empty_tseq_mlen_cnt_list(self):
        with self.assertRaises(ValueError):
            cml.fmt_tseq_rec_mlen_cnt(self.tseq_rec_list, [])

    def test_empty_both_input(self):
        with self.assertRaises(ValueError):
            cml.fmt_tseq_rec_mlen_cnt([], [])

    def test_diff_mlen(self):
        with self.assertRaises(AssertionError):
            cml.fmt_tseq_rec_mlen_cnt(self.tseq_rec_list, self.tseq_diff_mlen_cnt_list)

    def test_dup_tseq(self):
        with self.assertRaises(AssertionError):
            cml.fmt_tseq_rec_mlen_cnt(self.tseq_rec_list, 
                self.tseq_same_mlen_cnt_list + self.tseq_same_mlen_cnt_list[:2])

    def test_tseq_rec_seq_not_in_tseq(self):
        with self.assertRaises(ValueError):
            cml.fmt_tseq_rec_mlen_cnt(self.tseq_rec_non_exist_list, self.tseq_same_mlen_cnt_list)

        with self.assertRaises(ValueError):
            cml.fmt_tseq_rec_mlen_cnt(self.tseq_rec_non_exist_list + self.tseq_rec_list, 
                self.tseq_same_mlen_cnt_list)

    def test_output_format(self):
        o_tbl_str = cml.fmt_tseq_rec_mlen_cnt(self.tseq_rec_list, self.tseq_same_mlen_cnt_list)
        ref_tbl_str = '\
coord\tstrand\tseq\t0\t1\t2\t3\t4\t5\n\
0\t+\tX\t0\t0\t0\t0\t0\t0\n\
10\t-\tY\t0\t0\t0\t0\t0\t0\n\
20\t-\tATGtACG\t0\t0\t0\t0\t0\t0\n\
30\t+\tATTGtACG\t0\t0\t0\t0\t0\t0\n\
40\t+\tACCGtACG\t0\t0\t0\t0\t0\t0\n\
50\t-\tATGtACG\t0\t0\t0\t0\t0\t0\n\
60\t+\tACCGtACG\t0\t0\t0\t0\t0\t0\n\
70\t+\tACCGtACG\t0\t0\t0\t0\t0\t0\n'
        self.assertEqual(o_tbl_str, ref_tbl_str)


if __name__ == '__main__':
    unittest.main()

#!/usr/bin/env python2.7
import unittest
import env
import src.sequtil as su

class TestGetSeqMatchLen(unittest.TestCase):
    def all_matched(self):
        self.assertEqual(su.get_seq_match_len('ACGTT', 'ACGTT'), 5)
        self.assertEqual(su.get_seq_match_len('A', 'A'), 1)
        self.assertEqual(su.get_seq_match_len('AC', 'AC'), 2)
        self.assertEqual(su.get_seq_match_len('CAG', 'CAG'), 3)
        self.assertEqual(su.get_seq_match_len('cxaA', 'cxaA'), 4)

    def mismatched(self):
        self.assertEqual(su.get_seq_match_len('CCGTT', 'ACGTT'), 0)
        self.assertEqual(su.get_seq_match_len('A', 'C'), 0)
        self.assertEqual(su.get_seq_match_len('AC', 'AG'), 1)
        self.assertEqual(su.get_seq_match_len('CAG', 'CAa'), 2)
        self.assertEqual(su.get_seq_match_len('cxaA', 'cxax'), 3)

    def different_len(self):
        self.assertEqual(su.get_seq_match_len('A', 'AC'), 1)
        self.assertEqual(su.get_seq_match_len('A', 'ACG'), 1)
        self.assertEqual(su.get_seq_match_len('A', 'ACGT'), 1)
        self.assertEqual(su.get_seq_match_len('AC', 'ACGT'), 2)
        self.assertEqual(su.get_seq_match_len('ACGGG', 'ACGGGGGGGGG'), 5)

        self.assertEqual(su.get_seq_match_len('ACCGG', 'ACGGGGGGGGG'), 2)
        self.assertEqual(su.get_seq_match_len('AAGGG', 'ACGGGGGGGGG'), 1)
        self.assertEqual(su.get_seq_match_len('GACGGG', 'ACGGGGGGGGG'), 0)

    def empty_str(self):
        self.assertEqual(su.get_seq_match_len('', 'ACGGGGGGGGG'), 0)
        self.assertEqual(su.get_seq_match_len('', ''), 0)
        self.assertEqual(su.get_seq_match_len('', 'A'), 0)


class TestTargetSeqMatchLenCounter(unittest.TestCase):
    def test_counter(self):
        tseq_cnter = su.TargetSeqMatchLenCounter('ACGTACG')
        self.assertEqual(len(tseq_cnter._mlen_cnt_dict), 0)

        self.assertEqual(tseq_cnter._mlen_cnt_dict[0], 0)
        self.assertEqual(tseq_cnter._mlen_cnt_dict[0], 0)
        self.assertEqual(tseq_cnter._mlen_cnt_dict[1], 0)

        tseq_cnter.add_query_seq('CCGTACG')
        self.assertEqual(tseq_cnter._mlen_cnt_dict[0], 1)
        self.assertEqual(tseq_cnter._mlen_cnt_dict[1], 0)

        tseq_cnter.add_query_seq('ACGTACG')
        self.assertEqual(tseq_cnter._mlen_cnt_dict[0], 1)
        self.assertEqual(tseq_cnter._mlen_cnt_dict[1], 0)
        self.assertEqual(tseq_cnter._mlen_cnt_dict[7], 1)

        tseq_cnter.add_query_seq('')
        self.assertEqual(tseq_cnter._mlen_cnt_dict[0], 2)
        self.assertEqual(tseq_cnter._mlen_cnt_dict[1], 0)
        self.assertEqual(tseq_cnter._mlen_cnt_dict[7], 1)


        tseq_cnter.add_query_seq('ACATGACBDSP')
        self.assertEqual(tseq_cnter._mlen_cnt_dict[0], 2)
        self.assertEqual(tseq_cnter._mlen_cnt_dict[1], 0)
        self.assertEqual(tseq_cnter._mlen_cnt_dict[2], 1)
        self.assertEqual(tseq_cnter._mlen_cnt_dict[7], 1)

        tseq_cnter.add_query_seq('ACGTACG')
        tseq_cnter.add_query_seq('ACGTACG')
        tseq_cnter.add_query_seq('ACGTACG')
        tseq_cnter.add_query_seq('ACGTACG')
        tseq_cnter.add_query_seq('ACGTACG')
        self.assertEqual(tseq_cnter._mlen_cnt_dict[7], 6)
        self.assertEqual(tseq_cnter._mlen_cnt_dict[0], 2)
        self.assertEqual(tseq_cnter._mlen_cnt_dict[1], 0)
        self.assertEqual(tseq_cnter._mlen_cnt_dict[2], 1)


if __name__ == '__main__':
    unittest.main()
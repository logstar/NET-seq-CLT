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





if __name__ == '__main__':
    unittest.main()
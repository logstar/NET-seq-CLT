#!/usr/bin/env python2.7
import unittest
import env
import src.dnautil as du

class TestComplement(unittest.TestCase):
    def test_normal_dna_seq(self):
        self.assertEqual(du.complement('AACCGGTT'), 'TTGGCCAA')
        self.assertEqual(du.complement('A'), 'T')
        self.assertEqual(du.complement('T'), 'A')
        self.assertEqual(du.complement('C'), 'G')
        self.assertEqual(du.complement('G'), 'C')

    def test_abnormal_dna_seq(self):
        self.assertEqual(du.complement(''), '')
        self.assertEqual(du.complement('xACG'), 'xTGC')
        self.assertEqual(du.complement('\t'), '\t')

class TestReverseComplement(unittest.TestCase):
    def test_normal_dna_seq(self):
        self.assertEqual(du.reverse_complement('AACCGGTT'), 'AACCGGTT')
        self.assertEqual(du.reverse_complement('A'), 'T')
        self.assertEqual(du.reverse_complement('T'), 'A')
        self.assertEqual(du.reverse_complement('C'), 'G')
        self.assertEqual(du.reverse_complement('G'), 'C')

    def test_abnormal_dna_seq(self):
        self.assertEqual(du.reverse_complement(''), '')
        self.assertEqual(du.reverse_complement('xACG'), 'CGTx')
        self.assertEqual(du.reverse_complement('A\t'), '\tT')




if __name__ == '__main__':
    unittest.main()
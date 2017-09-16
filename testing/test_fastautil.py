#!/usr/bin/env python2.7
import unittest
import env
import src.fastautil as fa
import src.dnautil as dnautil

class TestSingleEntryFastaFile(unittest.TestCase):
    def test_seq_len(self):
        fa_file = fa.SingleEntryFastaFile('data/NC_000913_3.fna')
        self.assertEqual(fa_file.ref_seq_len, 4641652)

        tfa_fn = 'data/test_fautil1.fna'
        with open(tfa_fn, 'w') as tfa_file:
            tfa_file.write('>Test\n')
            tfa_file.write('A' * 50 + '\n')
            tfa_file.write('C' * 16 + '\n')
        test_fa_instance = fa.SingleEntryFastaFile(tfa_fn)
        self.assertEqual(test_fa_instance.ref_seq_len, 66)

    def test_get_seq(self):
        tfa_fn = 'data/test_fautil2.fna'
        with open(tfa_fn, 'w') as tfa_file:
            tfa_file.write('>Test\n')
            tfa_file.write('ACGTA' * 10 + '\n')
            tfa_file.write('C' * 16 + '\n')
        test_fa_instance = fa.SingleEntryFastaFile(tfa_fn)
        self.assertEqual(test_fa_instance.get_seq(0, 5), 'ACGTA')
        self.assertEqual(test_fa_instance.get_seq(5, 10), 'ACGTA')
        self.assertEqual(test_fa_instance.get_seq(3, 10), 'TAACGTA')
        self.assertEqual(test_fa_instance.get_seq(50, 55), 'CCCCC')
        self.assertEqual(test_fa_instance.get_seq(45, 55), 'ACGTACCCCC')

    def test_get_circularized_seq(self):
        tfa_fn = 'data/test_fautil2.fna'
        with open(tfa_fn, 'w') as tfa_file:
            tfa_file.write('>Test\n')
            tfa_file.write('ACGTA' * 10 + '\n')
            tfa_file.write('C' * 16 + '\n')
        test_fa_instance = fa.SingleEntryFastaFile(tfa_fn)

        self.assertEqual(test_fa_instance.get_circularized_seq(0, 5), 'ACGTA')
        self.assertEqual(test_fa_instance.get_circularized_seq(5, 5), 'ACGTA')
        self.assertEqual(test_fa_instance.get_circularized_seq(3, 7), 'TAACGTA')
        self.assertEqual(test_fa_instance.get_circularized_seq(50, 5), 'CCCCC')
        self.assertEqual(test_fa_instance.get_circularized_seq(45, 10), 'ACGTACCCCC')

        self.assertEqual(test_fa_instance.get_circularized_seq(45, 21), 'ACGTA' + 'C' * 16)
        self.assertEqual(test_fa_instance.get_circularized_seq(45, 22), 'ACGTA' + 'C' * 16 + 'A')
        self.assertEqual(test_fa_instance.get_circularized_seq(45, 23), 'ACGTA' + 'C' * 16 + 'AC')
        self.assertEqual(test_fa_instance.get_circularized_seq(45, 24), 'ACGTA' + 'C' * 16 + 'ACG')
        self.assertEqual(test_fa_instance.get_circularized_seq(45, 25), 'ACGTA' + 'C' * 16 + 'ACGT')
        self.assertEqual(test_fa_instance.get_circularized_seq(45, 26), 'ACGTA' + 'C' * 16 + 'ACGTA')
        self.assertEqual(test_fa_instance.get_circularized_seq(45, 21 + 66), 'ACGTA' + 'C' * 16 + test_fa_instance.ref_seq)

        self.assertEqual(test_fa_instance.get_circularized_seq(45, 21 + 67), 'ACGTA' + 'C' * 16 + test_fa_instance.ref_seq + 'A')
        self.assertEqual(test_fa_instance.get_circularized_seq(45, 21 + 68), 'ACGTA' + 'C' * 16 + test_fa_instance.ref_seq + 'AC')
        self.assertEqual(test_fa_instance.get_circularized_seq(45, 21 + 69), 'ACGTA' + 'C' * 16 + test_fa_instance.ref_seq + 'ACG')
        self.assertEqual(test_fa_instance.get_circularized_seq(45, 21 + 70), 'ACGTA' + 'C' * 16 + test_fa_instance.ref_seq + 'ACGT')
        self.assertEqual(test_fa_instance.get_circularized_seq(45, 21 + 71), 'ACGTA' + 'C' * 16 + test_fa_instance.ref_seq + 'ACGTA')

    def test_get_stranded_circularized_seq(self):
        tfa_fn = 'data/test_fautil2.fna'
        with open(tfa_fn, 'w') as tfa_file:
            tfa_file.write('>Test\n')
            tfa_file.write('ACGTA' * 10 + '\n')
            tfa_file.write('C' * 16 + '\n')
        test_fa_instance = fa.SingleEntryFastaFile(tfa_fn)
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(0, 0, 5 - 1, '+'), 'aCGTA')
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(5, 0, 5 - 1, '+'), 'aCGTA')
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(3, 0, 7 - 1, '+'), 'tAACGTA')
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(50, 0, 5 - 1, '+'), 'cCCCC')
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 0, 10 - 1, '+'), 'aCGTACCCCC')

        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 0, 21 - 1, '+'), 'aCGTA' + 'C' * 16)
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 0, 22 - 1, '+'), 'aCGTA' + 'C' * 16 + 'A')
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 0, 23 - 1, '+'), 'aCGTA' + 'C' * 16 + 'AC')
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 0, 24 - 1, '+'), 'aCGTA' + 'C' * 16 + 'ACG')
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 0, 25 - 1, '+'), 'aCGTA' + 'C' * 16 + 'ACGT')
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 0, 26 - 1, '+'), 'aCGTA' + 'C' * 16 + 'ACGTA')

        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 0, 21 + 66 - 1, '+'), 'aCGTA' + 'C' * 16 + test_fa_instance.ref_seq)
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 0, 21 + 67 - 1, '+'), 'aCGTA' + 'C' * 16 + test_fa_instance.ref_seq + 'A')
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 0, 21 + 68 - 1, '+'), 'aCGTA' + 'C' * 16 + test_fa_instance.ref_seq + 'AC')
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 0, 21 + 69 - 1, '+'), 'aCGTA' + 'C' * 16 + test_fa_instance.ref_seq + 'ACG')
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 0, 21 + 70 - 1, '+'), 'aCGTA' + 'C' * 16 + test_fa_instance.ref_seq + 'ACGT')
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 0, 21 + 71 - 1, '+'), 'aCGTA' + 'C' * 16 + test_fa_instance.ref_seq + 'ACGTA')

        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(0, 5 - 1, 0, '-'), dnautil.reverse_complement('aCGTA'))
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(5, 5 - 1, 0, '-'), dnautil.reverse_complement('aCGTA'))
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(3, 7 - 1, 0, '-'), dnautil.reverse_complement('tAACGTA'))
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(50, 5 - 1, 0, '-'), dnautil.reverse_complement('cCCCC'))
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 10 - 1, 0, '-'), dnautil.reverse_complement('aCGTACCCCC'))
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 21 - 1, 0, '-'), dnautil.reverse_complement('aCGTA' + 'C' * 16))
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 22 - 1, 0, '-'), dnautil.reverse_complement('aCGTA' + 'C' * 16 + 'A'))
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 23 - 1, 0, '-'), dnautil.reverse_complement('aCGTA' + 'C' * 16 + 'AC'))
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 24 - 1, 0, '-'), dnautil.reverse_complement('aCGTA' + 'C' * 16 + 'ACG'))
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 25 - 1, 0, '-'), dnautil.reverse_complement('aCGTA' + 'C' * 16 + 'ACGT'))
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 26 - 1, 0, '-'), dnautil.reverse_complement('aCGTA' + 'C' * 16 + 'ACGTA'))
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 21 + 66 - 1, 0, '-'), dnautil.reverse_complement('aCGTA' + 'C' * 16 + test_fa_instance.ref_seq))
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 21 + 67 - 1, 0, '-'), dnautil.reverse_complement('aCGTA' + 'C' * 16 + test_fa_instance.ref_seq + 'A'))
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 21 + 68 - 1, 0, '-'), dnautil.reverse_complement('aCGTA' + 'C' * 16 + test_fa_instance.ref_seq + 'AC'))
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 21 + 69 - 1, 0, '-'), dnautil.reverse_complement('aCGTA' + 'C' * 16 + test_fa_instance.ref_seq + 'ACG'))
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 21 + 70 - 1, 0, '-'), dnautil.reverse_complement('aCGTA' + 'C' * 16 + test_fa_instance.ref_seq + 'ACGT'))
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(45, 21 + 71 - 1, 0, '-'), dnautil.reverse_complement('aCGTA' + 'C' * 16 + test_fa_instance.ref_seq + 'ACGTA'))

        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(0, 4, 4, '+'), 'CCCCaCGTA')

        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(0, 4, 4, '-'), dnautil.reverse_complement('CCCCaCGTA'))
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(0, 4, 1, '-'), dnautil.reverse_complement('CaCGTA'))
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(0, 4, test_fa_instance.ref_seq_len, '-'), dnautil.reverse_complement(test_fa_instance.ref_seq + 'aCGTA'))
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(0, 4, 2 * test_fa_instance.ref_seq_len, '-'), dnautil.reverse_complement(2 * test_fa_instance.ref_seq + 'aCGTA'))

        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(0, 4, 1, '+'), 'CCCCaC')
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(0, 4, test_fa_instance.ref_seq_len, '+'), 'CCCCa' + test_fa_instance.ref_seq[1:] + 'A')
        self.assertEqual(test_fa_instance.get_stranded_circularized_seq(0, 4, 2 * test_fa_instance.ref_seq_len, '+'), 'CCCCa' + test_fa_instance.ref_seq[1:] + test_fa_instance.ref_seq + 'A')


if __name__ == '__main__':
    unittest.main()


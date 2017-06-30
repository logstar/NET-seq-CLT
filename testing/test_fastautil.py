#!/usr/bin/env python2.7
import unittest
import env
import src.fastautil as fa

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


if __name__ == '__main__':
    unittest.main()


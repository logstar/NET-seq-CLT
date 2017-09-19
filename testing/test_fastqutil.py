#!/usr/bin/env python2.7
import env
import unittest
from src import fastqutil

class TestFastqIterator(unittest.TestCase):
    def test_reading(self):
        seq_iter = fastqutil.iterate_fastq_file('data/linenum.txt')
        self.assertEqual(seq_iter.next(), ('1', '2', '3', '4'))
        self.assertEqual(seq_iter.next(), ('5', '6', '7', '8'))
        self.assertEqual(seq_iter.next(), ('9', '10', '11', '12'))
        with self.assertRaises(StopIteration):
            seq_iter.next()

        seq_iter = fastqutil.iterate_fastq_file('data/linenum.txt')
        self.assertEqual(seq_iter.next(), ('1', '2', '3', '4'))
        self.assertEqual(seq_iter.next(), ('5', '6', '7', '8'))
        self.assertEqual(seq_iter.next(), ('9', '10', '11', '12'))
        with self.assertRaises(StopIteration):
            seq_iter.next()

        seq_iter = fastqutil.iterate_fastq_file('data/empty_line.txt')
        self.assertEqual(seq_iter.next(), ('1', '', '3', ''))
        self.assertEqual(seq_iter.next(), ('5', '6', '7', '8'))
        with self.assertRaises(StopIteration):
            seq_iter.next()

        seq_iter = fastqutil.iterate_fastq_file('data/line_num_not_multiple_of_4.txt')
        self.assertEqual(seq_iter.next(), ('1', '', '3', ''))
        self.assertEqual(seq_iter.next(), ('5', '6', '7', '8'))
        with self.assertRaises(ValueError):
            seq_iter.next()


class TestIntToQScore(unittest.TestCase):
    def test_conversion(self):
        self.assertEqual(fastqutil.int_to_qscore(0), chr(33))
        self.assertEqual(fastqutil.int_to_qscore(1), chr(34))
        self.assertEqual(fastqutil.int_to_qscore(2), chr(35))
        self.assertEqual(fastqutil.int_to_qscore(3), chr(36))
        self.assertEqual(fastqutil.int_to_qscore(93), chr(33 + 93))

    def test_input_exceptions(self):
        with self.assertRaises(ValueError):
            fastqutil.int_to_qscore('1')

        with self.assertRaises(ValueError):
            fastqutil.int_to_qscore('a')

        with self.assertRaises(ValueError):
            fastqutil.int_to_qscore('!')

        with self.assertRaises(ValueError):
            fastqutil.int_to_qscore(-1)

        with self.assertRaises(ValueError):
            fastqutil.int_to_qscore(94)

        with self.assertRaises(ValueError):
            fastqutil.int_to_qscore(1.1)


if __name__ == '__main__':
    unittest.main()


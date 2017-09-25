#!/usr/bin/env python2.7
import unittest
import env
import src.trim_fastq as tf

class TestTrimFastq(unittest.TestCase):
    def test_result(self):
        i_fn = 'data/test_trim1.fastq'
        with open(i_fn, 'w') as i_file:
            # Good
            # 5' trim length 3
            # final length 7
            i_file.write(self.get_fastq_record('1', 'GGGAAAAAAAG', chr(50)))
            i_file.write(self.get_fastq_record('2', 'GGCCAAAAAAAC', chr(50)))
            i_file.write(self.get_fastq_record('3', 'GGGGAAAAAAGG', chr(50)))
            i_file.write(self.get_fastq_record('4', 'GGGTAAAAAACC', chr(50)))
            i_file.write(self.get_fastq_record('5', 'GGGCCAAAAATT', chr(50)))
            # short
            i_file.write(self.get_fastq_record('6', 'GGGCCAAAA', chr(50)))
            i_file.write(self.get_fastq_record('7', 'GGGCCAAA', chr(50)))
            # qual failed
            i_file.write(self.get_fastq_record('8', 'GGGCCAAAAATT', chr(49)))
            i_file.write(self.get_fastq_record('9', 'GGGCCAAAAATT', chr(49)))

        tf.trim_fastq([i_fn], 3, 7, chr(50), 'data/test_trim1_s.txt', 
                      'data/test_trim1_p.txt')
        with open('data/test_trim1_p.txt', 'r') as tmp_rfile:
            result_str = tmp_rfile.read()

        expected_result = ''.join((self.get_fastq_record('1', 'AAAAAAA', chr(50)),
                                   self.get_fastq_record('2', 'CAAAAAA', chr(50)),
                                   self.get_fastq_record('3', 'GAAAAAA', chr(50)),
                                   self.get_fastq_record('4', 'TAAAAAA', chr(50)),
                                   self.get_fastq_record('5', 'CCAAAAA', chr(50))))
        self.assertEqual(result_str, expected_result)

    @staticmethod
    def get_fastq_record(seqid, seq, qchar):
        qscore = qchar * len(seq)
        return "@%s\n%s\n+%s\n%s\n" % (seqid, seq, '', qscore)

if __name__ == '__main__':
    unittest.main()

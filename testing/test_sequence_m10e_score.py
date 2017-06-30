#!/usr/bin/env python2.7
import unittest
import env
import src.sequence_m10e_score as m10es


class TestSeqM10EScore(unittest.TestCase):
    def test_m10e_score_calculation(self):
        self.assertEqual(m10es.m10e_score('TATAAT'), 6)
        self.assertEqual(m10es.m10e_score('xAxxxT'), 2)
        self.assertEqual(m10es.m10e_score('TxTAAT'), 0)
        self.assertEqual(m10es.m10e_score('TATAAx'), 0)
        self.assertEqual(m10es.m10e_score('ACGTAG'), 0)
        self.assertEqual(m10es.m10e_score('TATACT'), 5)
        self.assertEqual(m10es.m10e_score('AATACT'), 4)
        self.assertEqual(m10es.m10e_score('AACACT'), 3)

    def test_incorrect_length(self):
        with self.assertRaises(ValueError):
            m10es.m10e_score('')
        with self.assertRaises(ValueError):
            m10es.m10e_score('A')
        with self.assertRaises(ValueError):
            m10es.m10e_score('AC')
        with self.assertRaises(ValueError):
            m10es.m10e_score('ACG')
        with self.assertRaises(ValueError):
            m10es.m10e_score('ACGT')
        with self.assertRaises(ValueError):
            m10es.m10e_score('ACGTA')
        with self.assertRaises(ValueError):
            m10es.m10e_score('ACGTACT')


class TestConsecutiveM10EScore(unittest.TestCase):
    def test_more_than_6bp(self):
        self.assertEqual(m10es.consecutive_m10e_score('TATAATT'), [6, 0])
        self.assertEqual(m10es.consecutive_m10e_score('xAxxxTATAAT'), [2, 0, 0, 0, 0, 6])
        self.assertEqual(m10es.consecutive_m10e_score('TxTAATAT'), [0, 0, 4])

        seq = 'A' * 20
        self.assertEqual(m10es.consecutive_m10e_score(seq), [0] * (len(seq) - 5))

    def test_6bp(self):
        self.assertEqual(m10es.consecutive_m10e_score('TATAAT'), [6])
        self.assertEqual(m10es.consecutive_m10e_score('xAxxxT'), [2])
        self.assertEqual(m10es.consecutive_m10e_score('TxTAAT'), [0])
        self.assertEqual(m10es.consecutive_m10e_score('TATAAx'), [0])
        self.assertEqual(m10es.consecutive_m10e_score('ACGTAG'), [0])
        self.assertEqual(m10es.consecutive_m10e_score('TATACT'), [5])
        self.assertEqual(m10es.consecutive_m10e_score('AATACT'), [4])
        self.assertEqual(m10es.consecutive_m10e_score('AACACT'), [3])

    def test_less_than_6bp(self):
        self.assertEqual(m10es.consecutive_m10e_score('TATAA'), [])
        self.assertEqual(m10es.consecutive_m10e_score('xxxx'), [])
        self.assertEqual(m10es.consecutive_m10e_score('TAC'), [])
        self.assertEqual(m10es.consecutive_m10e_score('AA'), [])
        self.assertEqual(m10es.consecutive_m10e_score('1'), [])
        self.assertEqual(m10es.consecutive_m10e_score(''), [])


class TestFrameFmtM10EScore(unittest.TestCase):
    def test_more_than_6bp(self):
        self.assertEqual(m10es.frame_fmt_m10e_score('TATAATT', 6), [(0, 'TATAAT', (6, )), 
                                                                    (1, 'ATAATT', (0, ))])
        self.assertEqual(m10es.frame_fmt_m10e_score('TATAATT', 7), [(0, 'TATAATT', (6, 0))])
        self.assertEqual(m10es.frame_fmt_m10e_score('TATAATT', 8), [])
        self.assertEqual(m10es.frame_fmt_m10e_score('TATAATT', 9), [])
        self.assertEqual(m10es.frame_fmt_m10e_score('TATAATT', 10), [])

        self.assertEqual(m10es.frame_fmt_m10e_score('xAxxxTATAAT', 6), [(0, 'xAxxxT', (2, )), 
                                                                        (1, 'AxxxTA', (0, )),
                                                                        (2, 'xxxTAT', (0, )),
                                                                        (3, 'xxTATA', (0, )),
                                                                        (4, 'xTATAA', (0, )),
                                                                        (5, 'TATAAT', (6, ))])

        self.assertEqual(m10es.frame_fmt_m10e_score('xAxxxTATAAT', 7), [(0, 'xAxxxTA', (2, 0)), 
                                                                        (1, 'AxxxTAT', (0, 0)),
                                                                        (2, 'xxxTATA', (0, 0)),
                                                                        (3, 'xxTATAA', (0, 0)),
                                                                        (4, 'xTATAAT', (0, 6))])

        self.assertEqual(m10es.frame_fmt_m10e_score('xAxxxTATAAT', 8), [(0, 'xAxxxTAT', (2, 0, 0)), 
                                                                        (1, 'AxxxTATA', (0, 0, 0)),
                                                                        (2, 'xxxTATAA', (0, 0, 0)),
                                                                        (3, 'xxTATAAT', (0, 0, 6))])

        seq = 'A' * 20
        frm_result = m10es.frame_fmt_m10e_score(seq, 6)
        self.assertEqual(map(lambda tup: tup[2], frm_result), [(0, )] * (len(seq) - 5))
        self.assertEqual(map(lambda tup: tup[1], frm_result), ['AAAAAA'] * (len(seq) - 5))

    def test_6bp(self):
        self.assertEqual(m10es.frame_fmt_m10e_score('TATAAT', 6), [(0, 'TATAAT', (6, ))])
        self.assertEqual(m10es.frame_fmt_m10e_score('xAxxxT', 6), [(0, 'xAxxxT', (2, ))])
        self.assertEqual(m10es.frame_fmt_m10e_score('TxTAAT', 6), [(0, 'TxTAAT', (0, ))])
        self.assertEqual(m10es.frame_fmt_m10e_score('TATAAx', 6), [(0, 'TATAAx', (0, ))])
        self.assertEqual(m10es.frame_fmt_m10e_score('ACGTAG', 6), [(0, 'ACGTAG', (0, ))])
        self.assertEqual(m10es.frame_fmt_m10e_score('TATACT', 6), [(0, 'TATACT', (5, ))])
        self.assertEqual(m10es.frame_fmt_m10e_score('AATACT', 6), [(0, 'AATACT', (4, ))])
        self.assertEqual(m10es.frame_fmt_m10e_score('AACACT', 6), [(0, 'AACACT', (3, ))])

    def test_less_than_6bp(self):
        self.assertEqual(m10es.frame_fmt_m10e_score('TATAA', 6), [])
        self.assertEqual(m10es.frame_fmt_m10e_score('xxxx', 6), [])
        self.assertEqual(m10es.frame_fmt_m10e_score('TAC', 6), [])
        self.assertEqual(m10es.frame_fmt_m10e_score('AA', 6), [])
        self.assertEqual(m10es.frame_fmt_m10e_score('1', 6), [])
        self.assertEqual(m10es.frame_fmt_m10e_score('', 6), [])

    def test_invalid_frame_size(self):
        with self.assertRaises(ValueError):
            m10es.frame_fmt_m10e_score('TATAA', 5)
        with self.assertRaises(ValueError):
            m10es.frame_fmt_m10e_score('TATAA', 4)
        with self.assertRaises(ValueError):
            m10es.frame_fmt_m10e_score('TATAA', 3)
        with self.assertRaises(ValueError):
            m10es.frame_fmt_m10e_score('TATAA', 2)
        with self.assertRaises(ValueError):
            m10es.frame_fmt_m10e_score('TATAA', 1)
        with self.assertRaises(ValueError):
            m10es.frame_fmt_m10e_score('TATAA', 0)
        with self.assertRaises(ValueError):
            m10es.frame_fmt_m10e_score('TATAA', -1)

    
if __name__ == '__main__':
    unittest.main()


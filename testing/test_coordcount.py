#!/usr/bin/env python2.7
import unittest
import env
import src.coordcount as cc

class TestCoordCountTable(unittest.TestCase):
    def setUp(self):
        ref_fn = 'data/test_coord_cnt_tbl.fna'
        with open(ref_fn, 'w') as ref_file:
            ref_file.write('>Test\n')
            ref_file.write('ACGTA' * 10 + '\n')
            ref_file.write('C' * 16 + '\n')

        ctbl_fn = 'data/coord_count_table1.txt'
        with open(ctbl_fn, 'w') as ctbl_file:
            ctbl_file.write('\t'.join(('0', '0', '9')) + '\n')
            ctbl_file.write('\t'.join(('1', '1', '8')) + '\n')
            ctbl_file.write('\t'.join(('2', '2', '7')) + '\n')
            ctbl_file.write('\t'.join(('3', '3', '6')) + '\n')
            ctbl_file.write('\t'.join(('4', '4', '5')) + '\n')
            ctbl_file.write('\t'.join(('5', '5', '4')) + '\n')
            ctbl_file.write('\t'.join(('6', '6', '3')) + '\n')
            ctbl_file.write('\t'.join(('7', '7', '2')) + '\n')
            ctbl_file.write('\t'.join(('8', '8', '1')) + '\n')
            ctbl_file.write('\t'.join(('9', '9', '0')) + '\n')
            ctbl_file.write('\t'.join(('10', '10', '20')) + '\n')
            ctbl_file.write('\t'.join(('11', '11', '22')) + '\n')
            ctbl_file.write('\t'.join(('65', '100', '200')) + '\n')
            ctbl_file.write('\t'.join(('64', '300', '400')) + '\n')
            ctbl_file.write('\t'.join(('63', '500', '600')) + '\n')

        self.test_table = cc.CoordCountTable(ctbl_fn, ref_fn)

    def test_init(self):
        self.assertEqual(self.test_table.ref_genome.ref_seq_len, 66)

        self.assertEqual(self.test_table.coord_cnt_dict[0].fwd_cnt, 0)
        self.assertEqual(self.test_table.coord_cnt_dict[0].rev_cnt, 9)

        self.assertEqual(self.test_table.coord_cnt_dict[6].fwd_cnt, 6)
        self.assertEqual(self.test_table.coord_cnt_dict[6].rev_cnt, 3)

        self.assertEqual(self.test_table.coord_cnt_dict[65].fwd_cnt, 100)
        self.assertEqual(self.test_table.coord_cnt_dict[65].rev_cnt, 200)

        self.assertEqual(self.test_table.coord_cnt_dict[63].fwd_cnt, 500)
        self.assertEqual(self.test_table.coord_cnt_dict[63].rev_cnt, 600)

        self.assertEqual(self.test_table.coord_cnt_dict[11].fwd_cnt, 11)
        self.assertEqual(self.test_table.coord_cnt_dict[11].rev_cnt, 22)

    def test_lr_range(self):
        self.assertEqual(self.test_table.get_left_right_range(10, '+', 5, 6), 
                         ([5, 6, 7, 8, 9], [11, 12, 13, 14, 15, 16]))
        self.assertEqual(self.test_table.get_left_right_range(10, '-', 5, 6), 
                         ([4, 5, 6, 7, 8, 9], [11, 12, 13, 14, 15]))

        self.assertEqual(self.test_table.get_left_right_range(2, '+', 5, 6), 
                         ([63, 64, 65, 0, 1], [3, 4, 5, 6, 7, 8]))
        self.assertEqual(self.test_table.get_left_right_range(63, '+', 5, 6), 
                         ([58, 59, 60, 61, 62], [64, 65, 0, 1, 2, 3]))

        self.assertEqual(self.test_table.get_left_right_range(2, '-', 6, 5), 
                         ([63, 64, 65, 0, 1], [3, 4, 5, 6, 7, 8]))
        self.assertEqual(self.test_table.get_left_right_range(63, '-', 6, 5), 
                         ([58, 59, 60, 61, 62], [64, 65, 0, 1, 2, 3]))

    def test_get_up_down_stream_seq_cnt_tuple(self):
        self.assertEqual(self.test_table.get_up_down_stream_seq_cnt_tuple(10, '+', 5, 6),
                         ('ACGTAaCGTAAC', (5, 6, 7, 8, 9, 10, 11, 0, 0, 0, 0, 0)))
        self.assertEqual(self.test_table.get_up_down_stream_seq_cnt_tuple(10, '-', 5, 6),
                         ('TTACGtTACGTT', (0, 0, 0, 0, 22, 20, 0, 1, 2, 3, 4, 5)))

        self.assertEqual(self.test_table.get_up_down_stream_seq_cnt_tuple(2, '+', 5, 6),
                         ('CCCACgTAACGT', (500, 300, 100, 0, 1, 2, 3, 4, 5, 6, 7, 8)))
        self.assertEqual(self.test_table.get_up_down_stream_seq_cnt_tuple(63, '+', 5, 6),
                         ('CCCCCcCCACGT', (0, 0, 0, 0, 0, 500, 300, 100, 0, 1, 2, 3)))
        # 'CCCCACgTAACG'
        self.assertEqual(self.test_table.get_up_down_stream_seq_cnt_tuple(2, '-', 5, 6),
                         ('CGTTAcGTGGGG', (2, 3, 4, 5, 6, 7, 8, 9, 200, 400, 600, 0)))
        # 'CCCCCCcCCACG'
        self.assertEqual(self.test_table.get_up_down_stream_seq_cnt_tuple(63, '-', 5, 6),
                         ('CGTGGgGGGGGG', (7, 8, 9, 200, 400, 600, 0, 0, 0, 0, 0, 0)))

    def test_rev_comp(self):
        self.assertEqual('GTTACGtTACGT', cc.CoordCountTable.rev_comp('ACGTAaCGTAAC'))
        self.assertEqual('A', cc.CoordCountTable.rev_comp('T'))
        self.assertEqual('ACGT', cc.CoordCountTable.rev_comp('ACGT'))
        self.assertEqual('TAGC', cc.CoordCountTable.rev_comp('GCTA'))
        self.assertEqual('accggt', cc.CoordCountTable.rev_comp('accggt'))


if __name__ == '__main__':
    unittest.main()

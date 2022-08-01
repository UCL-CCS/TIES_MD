#!/usr/bin/env python

import unittest
import numpy as np
from TIES_MD.ties_analysis.engines.namd import *


class Test_NAMD(unittest.TestCase):

    def test_alch_read(self):

        old_data = np.array([[ -12.4162,  -12.2928,   -8.2692,    2.5405,  -13.4956,  -12.4687],
                             [ 111.7324,   74.877,    72.9064,   71.3811,   80.7265,   77.0065],
                             [ -26.7703,  -63.6872,  -82.702,  -74.6786,  -87.6489,  -62.0604],
                             [-113.6954, -117.6372, -111.017,  -114.1721, -117.8366, -116.6829]])

        new_data = np.array([[ -11.5503,  -10.8766,    0.8469,   -8.9049,   22.1578,    3.1568],
                             [  99.995,    76.5626,   99.1309,   94.1263,   57.5984,   67.8794],
                             [ -70.0197,  -32.4018,  -67.8226,  -37.9828,  -87.3693,   -1.7703],
                             [-123.728,  -106.6641, -114.4813, -129.2575, -104.4866, -114.1935]])

        data_base = './data/test_read_data'

        old = read_alch_file(os.path.join(data_base, 'old.alch'), 2.11, 6)
        new = read_alch_file(os.path.join(data_base, 'new.alch'), 3, 6)

        self.assertIsNone(np.testing.assert_array_equal(old, old_data, 'Reading NAMD < 2.12 alch files incorrectly'))
        self.assertIsNone(np.testing.assert_array_equal(new, new_data, 'Reading NAMD > 2.12 alch files incorrectly'))

    def test_collating(self):
        vdw_a = []
        vdw_d = []
        ele_a = []
        ele_d = []

        namd_ti = NAMD('TI', None, None, None, None, [0], vdw_a, vdw_d, ele_a, ele_d, 2.12)
        data_dir = './data/test_collate_namd/l6-l14/lig'

        ti_ans = np.load(os.path.join(data_dir, 'assembled_ti.npy'))
        ti_built = namd_ti.collate_data('./data', 'test_collate_namd', 'l6-l14', 'lig')

        self.assertIsNone(np.testing.assert_array_equal(ti_ans, ti_built, 'TI result built incorrectly for NAMD2'))


if __name__ == '__main__':
    unittest.main()

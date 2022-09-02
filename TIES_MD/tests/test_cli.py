#!/usr/bin/env python

import unittest
from TIES_MD.cli import *


class Test_CLI(unittest.TestCase):

    def test_read_cfg(self):
        cfg_file = './systems/hydration/TIES.cfg'
        arg_dict = read_config(cfg_file)

        exp_dict = {'sampling_per_window': '2*unit.picosecond',
                    'equili_per_window': '2*unit.picoseconds',
                    'methods': 'FEP,TI',
                    'repeats': '2',
                    'elec_edges': '0.5,1.0',
                    'ster_edges': '0.0,0.5',
                    'windows': '3',
                    'constraint_file': 'na',
                    'constraint_column': 'beta_factor',
                    'input_type': 'AMBER',
                    'cell_basis_vec1': '41.5453,0.0,0.0',
                    'cell_basis_vec2': '0.0,45.082572999999995,0.0',
                    'cell_basis_vec3': '0.0,0.0,42.04022'}

        for k in exp_dict.keys():
            self.assertEqual(arg_dict[k], exp_dict[k], "Failed to read arg {} in config file".format(k))


if __name__ == '__main__':
    unittest.main()


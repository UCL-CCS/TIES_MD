#!/usr/bin/env python

import unittest
import shutil
from simtk import unit
from simtk.openmm import Vec3

from TIES_MD.TIES import *

class Test_TIES(unittest.TestCase):

    def test_initilization(self):
        args_dict = {'engine': 'openmm',
                     'temperature': '300*unit.kelvin',
                     'pressure': '1*unit.atmospheres',
                     'sampling_per_window': '2*unit.picosecond',
                     'equili_per_window': '2*unit.picoseconds',
                     'methods': 'FEP,TI',
                     'total_reps': '2',
                     'split_run': '0',
                     'elec_edges': '0.5,1.0',
                     'ster_edges': '0.0,0.5',
                     'global_lambdas': '0.0,0.5,1.0',
                     'constraint_file': 'cons.pdb',
                     'constraint_column': 'beta_factor',
                     'box_type': 'na',
                     'input_type': 'AMBER',
                     'cell_basis_vec1': '46.644591,0.0,0.0',
                     'cell_basis_vec2': '0.0,46.888166,0.0',
                     'cell_basis_vec3': '0.0,0.0,41.59781000000001'}

        test_msg = '{} was not initialized correctly.'

        test_init = TIES(cwd='./', run_type='class', exp_name='sol', devices=[0], node_id='test',
             windows_mask=[0,2], periodic=True, lam=None, **args_dict)

        self.assertEqual(test_init.sampling_per_window, unit.Quantity(2, unit.picoseconds),
                         test_msg.format('Sampling per window'))
        self.assertEqual(test_init.equili_per_window, unit.Quantity(2, unit.picoseconds),
                         test_msg.format('Equilibration per window'))

        self.assertEqual(test_init.elec_edges, [0.5, 1.0],
                         test_msg.format('Lambda schedule electric edges'))
        self.assertEqual(test_init.ster_edges, [0.0, 0.5],
                         test_msg.format('Lambda schedule steric edges'))

        self.assertEqual(test_init.methods, ['FEP', 'TI'],
                         test_msg.format('Methods'))
        self.assertEqual(test_init.total_reps, 2,
                         test_msg.format('Total number of repeats'))
        self.assertEqual(test_init.split_run, 0,
                         test_msg.format('Split run flag'))
        self.assertEqual(test_init.windows, 3,
                         test_msg.format('Number of windows'))

        self.assertEqual(test_init.constraint_file, 'cons.pdb',
                         test_msg.format('Constraint file location'))
        self.assertEqual(test_init.constraint_column, 'beta_factor',
                         test_msg.format('Constraint column'))

        self.assertEqual(test_init.basis_vectors, [Vec3(*[46.644591,0.0,0.0])*unit.angstrom,
                                                   Vec3(*[0.0,46.888166,0.0])*unit.angstrom,
                                                   Vec3(*[0.0,0.0,41.59781000000001])*unit.angstrom],
                         test_msg.format('Box size basis vectors'))

        self.assertEqual(test_init.absolute, False, test_msg.format('Absolute or relative option'))
        self.assertEqual(test_init.input_type, 'AMBER', test_msg.format('Input type'))
        self.assertEqual(test_init.windows_mask, [0,2], test_msg.format('Windows mask'))
        self.assertEqual(test_init.devices, [0], test_msg.format('Number of GPUs'))
        self.assertEqual(test_init.exp_name, 'sol', test_msg.format('Experiment name'))
        self.assertEqual(test_init.periodic, True, test_msg.format('Spacial periodicity flag'))

if __name__ == '__main__':
    unittest.main()

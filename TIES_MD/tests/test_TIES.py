#!/usr/bin/env python
import os
import unittest
import shutil
from simtk import unit
from simtk.openmm import Vec3
import io
from pathlib import Path
import glob
import copy
import site

from TIES_MD.TIES import *

test_args = {'engine': 'openmm',
             'temperature': '300.0*unit.kelvin',
             'pressure': '1.0*unit.atmospheres',
             'sampling_per_window': '0.002*unit.nanoseconds',
             'equili_per_window': '0.002*unit.nanoseconds',
             'methods': 'FEP,TI',
             'total_reps': '3',
             'split_run': '0',
             'elec_edges': '0.5,1.0',
             'ster_edges': '0.0,0.5',
             'global_lambdas': '0.0,0.5,1.0',
             'constraint_file': 'cons.pdb',
             'constraint_column': 'beta_factor',
             'input_type': 'AMBER',
             'cell_basis_vec1': '46.644591,0.0,0.0',
             'cell_basis_vec2': '0.0,46.888166,0.0',
             'cell_basis_vec3': '0.0,0.0,41.59781000000001'}

class Test_TIES(unittest.TestCase):

    def test_initialization(self):
        args_dict = copy.deepcopy(test_args)

        test_msg = '{} was not initialized correctly.'

        test_init = TIES(cwd='./', run_type='class', exp_name='sol', devices=[0], rep_id=None,
                         windows_mask=[0, 2], periodic=True, lam=None, **args_dict)

        self.assertEqual(test_init.sampling_per_window, unit.Quantity(0.002, unit.nanoseconds),
                         test_msg.format('Sampling per window'))
        self.assertEqual(test_init.equili_per_window, unit.Quantity(0.002, unit.nanoseconds),
                         test_msg.format('Equilibration per window'))

        self.assertEqual(test_init.elec_edges, [0.5, 1.0],
                         test_msg.format('Lambda schedule electric edges'))
        self.assertEqual(test_init.ster_edges, [0.0, 0.5],
                         test_msg.format('Lambda schedule steric edges'))

        self.assertEqual(test_init.methods, ['FEP', 'TI'],
                         test_msg.format('Methods'))
        self.assertEqual(test_init.total_reps, 3,
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

    #TODO add these tests for better coverage.
    def test_hpc_sub(self):
        pass

    def test_api(self):
        pass

    def test_analysis_cfg(self):
        pass

    def test_build_dirs(self):
        pass


    def test_namd_conf(self):

        #init a dummy input
        args_dict = copy.deepcopy(test_args)

        #make a output dir
        out_dir = os.path.join('data', 'test_replica-confs', 'tmp') #this is a relative such that CI dir structure is not
        #in generated conf scripts therfore we can match them.
        Path(out_dir).mkdir(parents=False, exist_ok=True)

        #iterate over the variables that effect how confs are written
        for engine in ['namd2.14', 'namd3']:
            for split, split_data in zip([0, 1], ['normal_run', 'split_run']):
                for cons, cons_data in zip(['cons.pdb', 'na'], ['cons', 'no_cons']):
                    args_dict['engine'] = engine
                    args_dict['split_run'] = split
                    args_dict['constraint_file'] = cons

                    test_init = TIES(cwd=out_dir, exp_name='sys_solv', devices=[0],
                                     rep_id=None, windows_mask=[0, 1], **args_dict)
                    #set up sim and write out some namd confs to check
                    test_init.setup()

                    #read base verified result
                    in_dir = os.path.join(os.getcwd(), 'data', 'test_replica-confs', engine, split_data, cons_data)
                    files_to_test = list(glob.iglob(os.path.join(in_dir, 'run*')))

                    for tst_path in files_to_test:
                        ref_path = tst_path.split(os.sep)[-1]
                        ref_path = os.path.join(out_dir, 'replica-confs', ref_path)

                        self.assertListEqual(list(io.open(tst_path)), list(io.open(ref_path)))

        #clean up out dir
        shutil.rmtree(out_dir)

    def test_update_cfg(self):

        #init a dummy input
        args_dict = copy.deepcopy(test_args)

        #make a output dir
        out_dir = os.path.join('data', 'tmp')
        Path(out_dir).mkdir(parents=False, exist_ok=True)

        # write TIES.cfg
        print(args_dict)
        test_init = TIES(cwd=out_dir, exp_name='sys_solv', devices=[0],
                         rep_id=None, windows_mask=[0, 1], **args_dict)
        test_init.update_cfg()

        #read cfg
        cfg = read_config(os.path.join(out_dir, 'TIES.cfg'))
        print(cfg)

        #test read cfg is same as input
        self.assertEqual(args_dict, cfg)

        # clean up out dir
        shutil.rmtree(out_dir)

if __name__ == '__main__':
    unittest.main()

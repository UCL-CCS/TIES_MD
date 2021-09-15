#!/usr/bin/env python

import unittest
import os
from TIES_MD.alch import *
from TIES_MD.lambdas import *
from simtk import unit
from TIES_MD.cli import read_config
from simtk.openmm import Vec3
import numpy as np
import shutil

from pathlib import Path

class Test_Alch(unittest.TestCase):

    def test_read_pdb(self):
        data_file_base = './systems/hydration/build'
        pdb = os.path.join(data_file_base, 'sol.pdb')
        lines = open(pdb).readlines()
        pdb_line = PDB_line(lines[5])

        fail = 'Failing to read pdb files'
        self.assertEqual(pdb_line.X, -4.250, fail)
        self.assertEqual(pdb_line.Y, -3.638, fail)
        self.assertEqual(pdb_line.Z, -3.031, fail)

        self.assertEqual(pdb_line.atom_serial_number, 4, fail)
        self.assertEqual(pdb_line.atom_name, ' H1 ', fail)
        self.assertEqual(pdb_line.residue_name, ' L1', fail)
        self.assertEqual(pdb_line.temperature_factor, -1.00, fail)
        self.assertEqual(pdb_line.occupancy, 1.00, fail)

    def test_get_alch_atoms(self):
        data_file_base = './systems/hydration/build'
        pdb = os.path.join(data_file_base, 'sol.pdb')
        app, dis = get_alchemical_atoms(pdb, absolute=False)
        fail = 'Misreading alchemical idxs in pdb file'
        self.assertEqual(dis, [2,3,4,5], fail)
        self.assertEqual(app, [32,33,34,35,36], fail)

    def test_get_constraints(self):
        data_file_base = './systems/hydration/build'
        pdb = os.path.join(data_file_base, 'cons.pdb')
        cons = get_constraints(pdb, 'occupancy')

        fail = 'Misreading constraint idxs in pdb file'

        self.assertEqual(cons, [unit.Quantity(1, unit.kilocalories_per_mole/unit.angstrom**2) for x in range(37)], fail)

    def test_grads(self):
        #define systems to test
        test_dirs = ['./systems/three_atom'
                     ]
        test_names = ['three'
                      ]

        config_files = ['TIES.cfg'
                        ]

        #define states to test
        states = [{'lambda_sterics_appear': 0, 'lambda_sterics_disappear': 1,
                   'lambda_electrostatics_appear': 0, 'lambda_electrostatics_disappear': 1},
                  {'lambda_sterics_appear': 0.5, 'lambda_sterics_disappear': 1.0,
                   'lambda_electrostatics_appear': 0.0, 'lambda_electrostatics_disappear': 0.5},
                  {'lambda_sterics_appear': 1.0, 'lambda_sterics_disappear': 0.5,
                   'lambda_electrostatics_appear': 0.5, 'lambda_electrostatics_disappear': 0.0},
                  {'lambda_sterics_appear': 1, 'lambda_sterics_disappear': 0,
                   'lambda_electrostatics_appear': 1, 'lambda_electrostatics_disappear': 0}]

        #Define know answers for each state, Calculated by hand from potentials
        known_vals = [{'lambda_sterics_appear': 0.263063638747951, 'lambda_sterics_disappear': 128.235461379091,
                       'lambda_electrostatics_appear': -30.9655122414654,
                       'lambda_electrostatics_disappear': -33.747881887578},
                      {'lambda_sterics_appear': 3.28224430120372, 'lambda_sterics_disappear': 128.2354613791,
                       'lambda_electrostatics_appear': -23.0239847106883,
                       'lambda_electrostatics_disappear': -24.4543486626991},
                      {'lambda_sterics_appear': 128.2354613791, 'lambda_sterics_disappear': 3.28224430120372,
                       'lambda_electrostatics_appear': -24.4240407261387,
                       'lambda_electrostatics_disappear': -23.0550629605553},
                      {'lambda_sterics_appear': 128.235461379091, 'lambda_sterics_disappear': 0.263063638747951,
                       'lambda_electrostatics_appear': -33.7765869123174,
                       'lambda_electrostatics_disappear': -31.0424648368013}]

        for cwd, exp_name, config in zip(test_dirs, test_names, config_files):


            #read arguments
            args_dict = read_config(os.path.join(cwd, config))

            temp = args_dict['temperature'].split('*unit.')
            temp = unit.Quantity(float(temp[0]), getattr(unit, temp[1]))

            pressure = args_dict['pressure'].split('*unit.')
            pressure = unit.Quantity(float(pressure[0]), getattr(unit, pressure[1]))

            cell_basis_vec1 = [float(x) for x in args_dict['cell_basis_vec1'].split(',')]
            cell_basis_vec2 = [float(x) for x in args_dict['cell_basis_vec2'].split(',')]
            cell_basis_vec3 = [float(x) for x in args_dict['cell_basis_vec3'].split(',')]

            #build system
            basis_vec = [Vec3(*cell_basis_vec1) * unit.angstrom,
                         Vec3(*cell_basis_vec2) * unit.angstrom,
                         Vec3(*cell_basis_vec3) * unit.angstrom]

            system = AlchSys(cwd, exp_name, temp, pressure, None, args_dict['constraint_column'], args_dict['methods'],
                             basis_vec, input_type='AMBER', absolute=False, periodic=True)
            NPT = system.build_simulation(system.NPT_alchemical_system, '0')
            NPT['sim'].context.setPositions(system.positions_file.positions)

            #Iterate over states and check the gradients are what we would expect and or agree with analytic calculations.
            for i, (param_vals_i, t_vals) in enumerate(zip(states, known_vals)):
                system.set_context_to_state(param_vals_i, NPT['sim'].context, NPT=True)

                num_g = {param: system.get_gradients(param, val, NPT['sim'].context, 0.0001).in_units_of(
                    unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole
                        for param, val in param_vals_i.items()}

                ana_g = {param: system.get_gradients(param, val, NPT['sim'].context, 0.0001, analitic_sterics=True).in_units_of(
                        unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole
                    for param, val in param_vals_i.items()}

                #Accuracy is bad here for analytic comparison
                self.assertEqual(round(num_g['lambda_sterics_appear'], 0), round(ana_g['lambda_sterics_appear'], 0),
                                 'Numerical and analyitic sterics do not agree for sys {}'.format(cwd))
                self.assertEqual(round(num_g['lambda_sterics_disappear'], 0), round(ana_g['lambda_sterics_disappear'], 0),
                                 'Numerical and analyitic sterics do not agree for sys {}'.format(cwd))

                #now compare to known values
                for k in param_vals_i.keys():
                    self.assertEqual(round(num_g[k], 6), round(t_vals[k], 6),
                                     '{} not being calculated correctly in state {}'.format(k, param_vals_i))

    def test_alch_system(self):
        # define systems to test
        test_dirs = ['./systems/hydration'
                     ]
        test_names = ['sol'
                      ]

        config_files = ['TIES.cfg'
                        ]

        # define states to test
        states = [{'lambda_sterics_appear': 1, 'lambda_sterics_disappear': 1,
               'lambda_electrostatics_appear': 1, 'lambda_electrostatics_disappear': 1},
              {'lambda_sterics_appear': 1, 'lambda_sterics_disappear': 0,
               'lambda_electrostatics_appear': 1, 'lambda_electrostatics_disappear': 0},
              {'lambda_sterics_appear': 0, 'lambda_sterics_disappear': 1,
               'lambda_electrostatics_appear': 0, 'lambda_electrostatics_disappear': 1},
              {'lambda_sterics_appear': 0, 'lambda_sterics_disappear': 0,
               'lambda_electrostatics_appear': 0, 'lambda_electrostatics_disappear': 0}]

        # Define know answers for each state, Calculated by hand from potentials
        known_vals = [-21694.272802303,
                      -21634.3811864154,
                      -21734.2116982551,
                      -21674.4343829331]

        for cwd, exp_name, config in zip(test_dirs, test_names, config_files):

            #read arguments
            args_dict = read_config(os.path.join(cwd, config))

            temp = args_dict['temperature'].split('*unit.')
            temp = unit.Quantity(float(temp[0]), getattr(unit, temp[1]))

            pressure = args_dict['pressure'].split('*unit.')
            pressure = unit.Quantity(float(pressure[0]), getattr(unit, pressure[1]))

            cell_basis_vec1 = [float(x) for x in args_dict['cell_basis_vec1'].split(',')]
            cell_basis_vec2 = [float(x) for x in args_dict['cell_basis_vec2'].split(',')]
            cell_basis_vec3 = [float(x) for x in args_dict['cell_basis_vec3'].split(',')]

            #build system
            basis_vec = [Vec3(*cell_basis_vec1) * unit.angstrom,
                         Vec3(*cell_basis_vec2) * unit.angstrom,
                         Vec3(*cell_basis_vec3) * unit.angstrom]

            system = AlchSys(cwd, exp_name, temp, pressure, None, args_dict['constraint_column'], args_dict['methods'],
                             basis_vec, input_type='AMBER', absolute=False, periodic=True)
            NPT = system.build_simulation(system.NPT_alchemical_system, '0')
            NPT['sim'].context.setPositions(system.positions_file.positions)

            #Iterate over states and check the energies are what we would expect
            for i, (param_vals_i, t_vals) in enumerate(zip(states, known_vals)):
                system.set_context_to_state(param_vals_i, NPT['sim'].context, NPT=True)
                energy  = NPT['sim'].context.getState(getEnergy=True).getPotentialEnergy().in_units_of(
                    unit.kilocalorie_per_mole)/unit.kilocalorie_per_mole
                self.assertEqual(round(energy, 2), round(t_vals, 2),
                                 'energy not being calculated correctly in state {}'.format(param_vals_i))


    def test_run_sim(self):
        # define systems to test
        test_dirs = ['./systems/three_atom'
                     ]
        test_names = ['three'
                      ]
        config_files = ['TIES.cfg'
                        ]

        for cwd, exp_name, config in zip(test_dirs, test_names, config_files):
            # read arguments
            args_dict = read_config(os.path.join(cwd, config))

            temp = args_dict['temperature'].split('*unit.')
            temp = unit.Quantity(float(temp[0]), getattr(unit, temp[1]))

            pressure = args_dict['pressure'].split('*unit.')
            pressure = unit.Quantity(float(pressure[0]), getattr(unit, pressure[1]))

            cell_basis_vec1 = [float(x) for x in args_dict['cell_basis_vec1'].split(',')]
            cell_basis_vec2 = [float(x) for x in args_dict['cell_basis_vec2'].split(',')]
            cell_basis_vec3 = [float(x) for x in args_dict['cell_basis_vec3'].split(',')]

            # build system
            basis_vec = [Vec3(*cell_basis_vec1) * unit.angstrom,
                         Vec3(*cell_basis_vec2) * unit.angstrom,
                         Vec3(*cell_basis_vec3) * unit.angstrom]

            system = AlchSys(cwd, exp_name, temp, pressure, None, args_dict['constraint_column'], args_dict['methods'],
                             basis_vec, input_type='AMBER', absolute=False, periodic=True)

            node_id = 0
            ids= System_ID(device_id='0', node_id=node_id)

            Lam = Lambdas([0.5, 1], [0.0, 0.5], [x/5 for x in range(0, 6)])
            mask= [0, 2]
            niter = 2
            equili_steps = 10
            steps_per_iter = 1

            #build some output dirs
            for d in ['results', 'simulation', 'equilibration']:
                for l in ['0', '1']:
                    path = os.path.join(cwd, 'LAMBDA_{}/rep{}/{}'.format(l, node_id, d))
                    Path(path).mkdir(parents=True, exist_ok=True)

            simulate_system(ids, system, Lam, mask, cwd, niter, equili_steps, steps_per_iter)

            #4 is number of lambda dimensions
            expected_shapes = [(4, niter), (len(Lam.schedule), niter)]

            #look for output
            for l in ['0', '1']:
                for ans, method in zip(expected_shapes, ['TI', 'FEP']):
                    path = os.path.join(cwd, 'LAMBDA_{}/rep{}/results/{}.npy'.format(l, node_id, method))
                    array = np.load(path)
                    self.assertEqual(array.shape, ans, 'Failed to generate results of expected shape')

            #delete output
            # build some output dirs
            for l in ['0', '1']:
                path = os.path.join(cwd, 'LAMBDA_{}'.format(l))
                shutil.rmtree(path)

if __name__ == '__main__':

    unittest.main()

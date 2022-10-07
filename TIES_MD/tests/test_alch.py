#!/usr/bin/env python
import copy
import unittest
import os
from TIES_MD.alch import *
from TIES_MD.lambdas import *
from TIES_MD.cli import read_config
import numpy as np
import shutil
from pathlib import Path

try:
    import openmm as mm
    from openmm import unit, app, Vec3
except ImportError:  # OpenMM < 7.6
    import simtk.openmm as mm
    from simtk.openmm import app
    from simtk import unit
    from simtk.openmm import Vec3

GLOBAL_PALT = 'CPU'

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
                             basis_vec, input_type='AMBER', absolute=False, periodic=True, platform=GLOBAL_PALT)
            NPT = system.build_simulation(system.NPT_alchemical_system, '0')
            NPT['sim'].context.setPositions(system.positions_file.positions)

            #Iterate over states and check the gradients agree with analytic calculations.
            for i, param_vals_i in enumerate(states):
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


    def test_alch_system(self):
        # define systems to test
        test_dirs = ['./systems/hydration'
                     ]
        test_names = ['sol'
                      ]

        config_files = ['TIES.cfg'
                        ]
        #outer list for systems inner lists for [appearing, disappearing]
        alchemical_atoms = [[[32, 33, 34, 35, 36], [2, 3, 4, 5]]
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

        for cwd, exp_name, config, atoms in zip(test_dirs, test_names, config_files, alchemical_atoms):
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

            ties_system = AlchSys(cwd, exp_name, temp, pressure, None, args_dict['constraint_column'], args_dict['methods'],
                             basis_vec, input_type='AMBER', absolute=False, periodic=True, platform=GLOBAL_PALT, debug=True)
            NPT = ties_system.build_simulation(ties_system.NPT_alchemical_system, '0')
            NPT['sim'].context.setPositions(ties_system .positions_file.positions)

            #create normal openmm system
            coord_file = os.path.join(cwd, 'build', exp_name + '.pdb')
            top_file = os.path.join(cwd, 'build', exp_name + '.prmtop')
            positions_file = mm.app.pdbfile.PDBFile(coord_file)
            topology_file = mm.app.AmberPrmtopFile(top_file)
            nonbondedMethod = app.PME
            openmm_system = topology_file.createSystem(nonbondedMethod=nonbondedMethod,
                                                 nonbondedCutoff=1.2 * unit.nanometers,
                                                 constraints=app.HBonds, rigidWater=True,
                                                 ewaldErrorTolerance=0.00001)
            openmm_system.setDefaultPeriodicBoxVectors(*basis_vec)
            openmm_system.addForce(mm.MonteCarloBarostat(pressure, temp, 25))
            null_cross_region_interactions(openmm_system, atoms[0], atoms[1])
            reduce_to_nonbonded(openmm_system)

            # grab nonbonded force from system
            for force_index, force in enumerate(openmm_system.getForces()):
                if isinstance(force, mm.NonbondedForce):
                    nonbonded_force = force

            # add switching function to nonbonded
            nonbonded_force.setSwitchingDistance(1.0 * unit.nanometers)
            nonbonded_force.setUseSwitchingFunction(True)
            nonbonded_force.setUseDispersionCorrection(False)

            #Iterate over states and check the energies are what we would expect
            for i, param_vals_i in enumerate(states):
                #use TIES
                ties_system.set_context_to_state(param_vals_i, NPT['sim'].context, NPT=True)
                energy_1 = NPT['sim'].context.getState(getEnergy=True).getPotentialEnergy().in_units_of(
                    unit.kilocalorie_per_mole)/unit.kilocalorie_per_mole

                #use OpenMM
                test_sys = copy.deepcopy(openmm_system)
                turn_off_interactions(param_vals_i, test_sys, atoms[0], atoms[1])
                friction = 0.3 / unit.picosecond
                time_step = 2.0 * unit.femtosecond
                integrator = mm.LangevinIntegrator(temp, friction, time_step)
                integrator.setConstraintTolerance(0.00001)
                platform = mm.Platform.getPlatformByName(GLOBAL_PALT)
                properties = {}
                sim = app.Simulation(topology_file.topology, test_sys, integrator, platform, properties)
                sim.context.setPositions(positions_file.positions)
                energy_2 = sim.context.getState(getEnergy=True).getPotentialEnergy().in_units_of(
                    unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole


                #check equal
                self.assertEqual(round(energy_1, 3), round(energy_2, 3),
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
                             basis_vec, input_type='AMBER', absolute=False, periodic=True, platform=GLOBAL_PALT)

            rep_id = 0
            ids = System_ID(device_id='0', rep_id=rep_id)

            Lam = Lambdas([0.5, 1], [0.0, 0.5], [x/5 for x in range(0, 6)])
            mask = [0, 2]
            niter = 2
            equili_steps = 10
            steps_per_iter = 1

            #build some output dirs
            for d in ['results', 'simulation', 'equilibration']:
                for l in ['0.00', '0.20']:
                    path = os.path.join(cwd, 'LAMBDA_{}/rep{}/{}'.format(l, rep_id, d))
                    Path(path).mkdir(parents=True, exist_ok=True)

            simulate_system(ids, system, Lam, mask, cwd, niter, equili_steps, steps_per_iter)

            #4 is number of lambda dimensions
            expected_shapes = [(4, niter), (len(Lam.schedule), niter)]

            #look for output
            for l in ['0.00', '0.20']:
                for ans, method in zip(expected_shapes, ['TI', 'FEP']):
                    path = os.path.join(cwd, 'LAMBDA_{}/rep{}/results/{}.npy'.format(l, rep_id, method))
                    array = np.load(path)
                    self.assertEqual(array.shape, ans, 'Failed to generate results of expected shape')

            #delete output
            for l in ['0.00', '0.20']:
                path = os.path.join(cwd, 'LAMBDA_{}'.format(l))
                shutil.rmtree(path)


def reduce_to_nonbonded(system):
    to_remove = []
    for force_index, force in enumerate(system.getForces()):
        if isinstance(force, mm.NonbondedForce):
            pass
        elif isinstance(force, mm.MonteCarloBarostat):
            pass
        else:
            to_remove.append(force_index)
    to_remove.sort()
    to_remove.reverse()
    for force_index in to_remove:
        system.removeForce(force_index)


def null_cross_region_interactions(system, appearing, disappearing):
    for force_index, force in enumerate(system.getForces()):
        if isinstance(force, mm.NonbondedForce):
            for atom1 in appearing:
                for atom2 in disappearing:
                    force.addException(atom1, atom2, 0.0, 1.0, 0.0, True)


def turn_off_interactions(state, system, appearing, disappearing):
    for force_index, force in enumerate(system.getForces()):
        if isinstance(force, mm.NonbondedForce):
            if state['lambda_sterics_appear'] == 0:
                for idx in range(force.getNumParticles()):
                    if idx in appearing:
                        q, sig, eps = force.getParticleParameters(idx)
                        force.setParticleParameters(idx, q, 1, 0)
                for idx in range(force.getNumExceptions()):
                    a, b, q, sig, eps = force.getExceptionParameters(idx)
                    if set((a, b)).intersection(appearing):
                        force.setExceptionParameters(idx, a, b, q, 1, 0)

            if state['lambda_sterics_disappear'] == 0:
                for idx in range(force.getNumParticles()):
                    if idx in disappearing:
                        q, sig, eps = force.getParticleParameters(idx)
                        force.setParticleParameters(idx, q, 1, 0)
                for idx in range(force.getNumExceptions()):
                    a, b, q, sig, eps = force.getExceptionParameters(idx)
                    if set((a, b)).intersection(disappearing):
                        force.setExceptionParameters(idx, a, b, q, 1, 0)

            if state['lambda_electrostatics_appear'] == 0:
                for idx in range(force.getNumParticles()):
                    if idx in appearing:
                        q, sig, eps = force.getParticleParameters(idx)
                        force.setParticleParameters(idx, 0, sig, eps)
                for idx in range(force.getNumExceptions()):
                    a, b, q, sig, eps = force.getExceptionParameters(idx)
                    if set((a, b)).intersection(appearing):
                        force.setExceptionParameters(idx, a, b, 0, sig, eps)

            if state['lambda_electrostatics_disappear'] == 0:
                for idx in range(force.getNumParticles()):
                    if idx in disappearing:
                        q, sig, eps = force.getParticleParameters(idx)
                        force.setParticleParameters(idx, 0, sig, eps)
                for idx in range(force.getNumExceptions()):
                    a, b, q, sig, eps = force.getExceptionParameters(idx)
                    if set((a, b)).intersection(disappearing):
                        force.setExceptionParameters(idx, a, b, 0, sig, eps)

if __name__ == '__main__':
    unittest.main()

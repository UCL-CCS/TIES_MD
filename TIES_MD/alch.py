#!/usr/bin/env python

__copyright__ = """
    Copyright 2021 Alexander David Wade
    This file is part of TIES MD

    TIES MD is free software: you can redistribute it and/or modify
    it under the terms of the Lesser GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    TIES MD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    Lesser GNU General Public License for more details.
    You should have received a copy of the Lesser GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

__license__ = "LGPL"

try:
    import openmm as mm
    from openmm import unit, app, Vec3
except ImportError:  # OpenMM < 7.6
    import simtk.openmm as mm
    from simtk.openmm import app
    from simtk import unit
    from simtk.openmm import Vec3

import numpy as np
import openmmtools

import os
import copy
import glob
import time
from sys import stdout

from .openmmtools.alchemy import ModifiedAbsoluteAlchemicalFactory, ModifiedAlchemicalState


class AlchSys(object):
    '''
    Class containing one leg of binding free energy calculation

    :param cwd: str, for current working directory
    :param name: str, name of system
    :param temperature: temperature of the simulation and target for the thermostat
    :param pressure: pressure of the simulation and target for the barostat
    :param constraint_file: str, for location of pdb which details constraints
    :param constraint_column: str, defines if constraints are read from the beta or occupancy columns
    :param methods: list of str, list for what methods we want to use example ['MBAR', 'TI']
    :param basis_vectors: list of list containing explicit cell vectors
    :param absolute: boolean determines if we are doing an absolute calculation
    :param periodic: boolean determines if we are doing PME or CutoffNonPeriodic
    :param platform: sting determines what platform OpenMM will target allowed values are ['CPU', 'CUDA', 'OpenCL']
    :param debug: bool, removes all forces but nonbonded for testing

    Note: GROMACS input and Absolute calculations are currently not tested.

    '''
    def __init__(self, cwd, name, temperature, pressure, constraint_file, constraint_column, methods,
                 basis_vectors, input_type='AMBER', absolute=False, periodic=True, platform='CUDA', debug=False):

        self.absolute = absolute
        self.platform = platform
        self.temp = temperature
        self.pressure = pressure

        ##System constants, these can be passed as arguments later
        self.friction = 0.3 / unit.picosecond
        self.time_step = 2.0 * unit.femtosecond

        if not os.path.exists(os.path.join(cwd, 'build')):
            raise FileNotFoundError('build directory {} found not to exist'.format(os.path.join(cwd, 'build')))

        self.methods = methods
        if constraint_file is not None:
            self.constraints = get_constraints(os.path.join(cwd, 'build', constraint_file), constraint_column)
        else:
            self.constraints = None

        if input_type == 'AMBER':
            # Using pdb for now could use inpcrd
            coord_file = os.path.join(cwd, 'build', name + '.pdb')
            alch_file = os.path.join(cwd, 'build', name + '.pdb')
            top_file = os.path.join(cwd, 'build', name + '.prmtop')
            #The positions in this file are used unless some alchemical atoms are fatally overlapped
            # see the function shift_alchemical_positions()
            self.positions_file = mm.app.pdbfile.PDBFile(coord_file)
            self.og_positions = self.positions_file.positions
            self.topology_file = mm.app.AmberPrmtopFile(top_file)

        elif input_type == 'GROMACS':
            # needs to be general
            gmx_dir = os.path.join(cwd, 'build')
            coord_file = os.path.join(cwd, 'build', name + '.gro')
            top_file = os.path.join(cwd, 'build', name + '.top')
            alch_file = os.path.join(cwd, 'build', name + '.pdb')
            # The positions in this file are used unless some alchemical atoms are fatally overlapped
            # see the function shift_alchemical_positions()
            self.positions_file = mm.app.GromacsGroFile(coord_file)
            self.og_positions = self.positions_file.positions
            #might be a problem for HPC where gromacs is not avaliable
            self.topology_file = mm.app.GromacsTopFile(top_file, periodicBoxVectors=self.positions_file.getPeriodicBoxVectors(),
                                                       includeDir=gmx_dir)
        else:
            raise TypeError('Unknown input type {} please chose from AMBER/GROMACS'.format(input_type))

        #Set perioic boundary conditions
        if periodic:
            nonbondedMethod = app.PME
        else:
            nonbondedMethod = app.CutoffNonPeriodic

        system = self.topology_file.createSystem(nonbondedMethod=nonbondedMethod , nonbondedCutoff=1.2 * unit.nanometers,
                                                 constraints=app.HBonds, rigidWater=True, ewaldErrorTolerance=0.00001)

        #grab nonbonded force from system
        for force_index, force in enumerate(system.getForces()):
            if isinstance(force, mm.NonbondedForce):
                nonbonded_force = force

        # add switching function to nonbonded
        nonbonded_force.setSwitchingDistance(1.0 * unit.nanometers)
        nonbonded_force.setUseSwitchingFunction(True)

        #Turn off dispertion correction to match namd
        nonbonded_force.setUseDispersionCorrection(False)

        #Extract indexs of appearing and disappearing atoms from input PDB file
        #If using gromacs how can we get this information?
        appear_idxs, disappear_idxs = get_alchemical_atoms(alch_file, self.absolute)

        if len(appear_idxs) == 0 or len(disappear_idxs) == 0:
            if self.absolute is False:
                raise ValueError('Performing relative calculation but only found one alchemical region in pdb')
            #Move alchemical atoms to dissappear region so we know where to expect them
            disappear_idxs = disappear_idxs+appear_idxs
            appear_idxs = []

        print('Appearing atoms {}'.format(appear_idxs))
        print('Disappearing atoms {}'.format(disappear_idxs))

        self.appearing_idxs = appear_idxs

        print('Default box vectors were:')
        print(system.getDefaultPeriodicBoxVectors())

        print('Switching to explicit box vectors:')
        system.setDefaultPeriodicBoxVectors(*basis_vectors)
        self.PBV = system.getDefaultPeriodicBoxVectors()
        for axis in self.PBV:
            print(axis)

        if debug:
            AlchSys.debug_force(self, system)

        # find angles, bonds and torsions which straddle alchemical regions, so we can turn them off later
        intersect_angles = AlchSys.get_intersect_angles(self, system, appear_idxs, disappear_idxs)
        intersect_bonds = AlchSys.get_intersect_bonds(self, system, appear_idxs, disappear_idxs)
        intersect_torsions = AlchSys.get_intersect_torsions(self, system, appear_idxs, disappear_idxs)

        self.rebuild_torsion(system, intersect_torsions)

        # create and alchemical system using OpenMMTools
        factory = ModifiedAbsoluteAlchemicalFactory(alchemical_pme_treatment='exact')

        disappear = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms=disappear_idxs, name='disappear',
                                                         annihilate_electrostatics=True, annihilate_sterics=True,
                                                         softcore_alpha=0.5, softcore_a=1, softcore_b=1, softcore_c=6)

        if not self.absolute:
            appear = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms=appear_idxs, alchemical_bonds=intersect_bonds,
                                                          alchemical_angles=intersect_angles,
                                                          name='appear', annihilate_electrostatics=True, annihilate_sterics=True,
                                                          softcore_alpha=0.5, softcore_a=1, softcore_b=1, softcore_c=6)

            # Create alchemical system. Note: by default appear and disappear cant see each other
            self.NVT_alchemical_system = factory.create_alchemical_system(system, alchemical_regions=[disappear, appear])

        else:
            # Create alchemical system with one region
            self.NVT_alchemical_system = factory.create_alchemical_system(system,
                                                                          alchemical_regions=[disappear])

        #Delcalare what derivatives to calculate
        for force in self.NVT_alchemical_system.getForces():
            if isinstance(force, mm.CustomNonbondedForce) or isinstance(force, mm.CustomBondForce):
                for i in range(force.getNumGlobalParameters()):
                    param_name = force.getGlobalParameterName(i)
                    if 'lambda_sterics' in param_name:
                        print('Adding derivative for {} to a CustomNonbondedForce or CustomBondForce'.format(param_name))
                        force.addEnergyParameterDerivative(param_name)
                        continue

        alchemical_state_d = ModifiedAlchemicalState.from_system(self.NVT_alchemical_system,
                                                                             parameters_name_suffix='disappear')

        if not self.absolute:
            alchemical_state_a = ModifiedAlchemicalState.from_system(self.NVT_alchemical_system,
                                                                             parameters_name_suffix='appear')
            composable_states = [alchemical_state_d, alchemical_state_a]
        else:
            composable_states = [alchemical_state_d]

        # compound alchemical states a and d into one state
        ts = openmmtools.states.ThermodynamicState(self.NVT_alchemical_system, self.temp)
        self.NVT_compound_state = openmmtools.states.CompoundThermodynamicState(thermodynamic_state=ts,
                                                                            composable_states=composable_states)
        #Turn off intersect angle, bond, torsions
        if len(intersect_angles) > 0:
            print('Found angles {} straddling alchemical regions these will be turned off'.format(intersect_angles))
            self.NVT_compound_state.lambda_angles_appear = 0.0
        if len(intersect_bonds) > 0:
            print('Found bonds {} straddling alchemical regions these will be turned off'.format(intersect_bonds))
            self.NVT_compound_state.lambda_bonds_appear = 0.0
        if len(intersect_torsions) > 0:
            print('Found torsion {} straddling alchemical regions these will be removed'.format(intersect_torsions))

        self.NPT_alchemical_system = copy.deepcopy(self.NVT_alchemical_system)
        self.NPT_compound_state = copy.deepcopy(self.NVT_compound_state)

        #if not periodic NPT is just a copy of NVT this is confusing as code but easy to implement for now.
        if periodic:
            #add barostat and set presssure
            self.NPT_alchemical_system.addForce(mm.MonteCarloBarostat(self.pressure, self.temp, 25))
            # Prime OpenMMtools to anticipate systems with barostats
            self.NPT_compound_state.pressure = self.pressure

        # Add constraints to NVT system
        if self.constraints is not None:
            print('Constraining NVT system')
            AlchSys.add_consraints(self)

    def amend_original_positions(self, positions):
        '''
        Function to update the stored positions if a clash of atoms is found during initialization.

        :param positions: list, for the positions of all atoms in the system to be used to initialise simulations
        '''
        self.og_positions = positions

    def test_sim(self, simulation):
        '''
        Function to test if system can be minimized, useful to find clashed atoms or undefined potentials

        :param simulation: OpenMM simulation object to test.
        '''
        simulation['sim'].context.setPositions(self.og_positions)
        mm.LocalEnergyMinimizer.minimize(simulation['sim'].context, maxIterations=10)


    def shift_alchemical_positions(self):
        '''
        Function to add a small pertubation to the positions of all appearing atoms in alchemical region.
        This can help to resolve any nans caused by overlapping atoms.
        '''

        new_pos = self.og_positions

        # shift atoms
        for i, pos in enumerate(new_pos):
            if i in self.appearing_idxs:
                new_pos[i] = pos + Vec3(0.00001, 0.00001, 0.00001)*unit.nanometers
                print('shifting atom {} from {} to {}'.format(i, pos, new_pos[i]))

        AlchSys.amend_original_positions(self, new_pos)

    def get_intersect_angles(self, system, appear_idxs, disappear_idxs):
        '''
        Function to get the idxs of angles which straddle the alchemical region.

        :param system: OpenMM system object
        :param appear_idxs: list of indexes for appearing atoms
        :param disappear_idxs: list of indexes for dissapearing atoms

        :return: list of indexes for angles which straddle alchemical regions
        '''
        intersect_angles = []
        if len(appear_idxs) == 0 or len(disappear_idxs) == 0:
            return intersect_angles
        for force_index, force in enumerate(system.getForces()):
            if isinstance(force, mm.HarmonicAngleForce):
                for idx in range(force.getNumAngles()):
                    a, b, c, theta, k = force.getAngleParameters(idx)
                    if set((a, b, c)).intersection(appear_idxs) and set((a, b, c)).intersection(disappear_idxs):
                        intersect_angles.append(idx)
        return intersect_angles

    def get_intersect_bonds(self, system, appear_idxs, disappear_idxs):
        '''
        Function to get the idxs of bonds which straddle the alchemical region.

        :param system: OpenMM system object
        :param appear_idxs: list of indexes for appearing atoms
        :param disappear_idxs: list of indexes for dissapearing atoms

        :return: list of indexes for bonds which straddle alchemical regions
        '''
        intersect_bonds = []
        if len(appear_idxs) == 0 or len(disappear_idxs) == 0:
            return intersect_bonds
        for force_index, force in enumerate(system.getForces()):
            if isinstance(force, mm.HarmonicBondForce):
                for idx in range(force.getNumBonds()):
                    a, b, r0, k = force.getBondParameters(idx)
                    if set((a, b)).intersection(appear_idxs) and set((a, b)).intersection(disappear_idxs):
                        intersect_bonds.append(idx)
        return intersect_bonds

    def get_intersect_torsions(self, system, appear_idxs, disappear_idxs):
        '''
        Function to get the idxs of torsions which straddle the alchemical region.

        :param system: OpenMM system object
        :param appear_idxs: list of indexes for appearing atoms
        :param disappear_idxs: list of indexes for disappearing atoms
        :return: list of indexes for torsions which straddle alchemical regions
        '''
        intersect_tor= []
        if len(appear_idxs) == 0 or len(disappear_idxs) == 0:
            return intersect_tor
        for force_index, force in enumerate(system.getForces()):
            if isinstance(force, mm.PeriodicTorsionForce):
                for idx in range(force.getNumTorsions()):
                    a, b, c, d, period, phase, k = force.getTorsionParameters(idx)
                    if set((a, b, c, d)).intersection(appear_idxs) and set((a, b, c, d)).intersection(disappear_idxs):
                        intersect_tor.append(idx)
        return intersect_tor

    def rebuild_torsion(self, system, intersect_torsions):
        '''
        Function to rebuild torsion force without any non physical torsion which straddle the alchemical region. We
        want to remove these fully as some times they result in nan evals.

        :param system: OpenMM system
        :param intersect_torsions: List of ints, ints references torsions that straddle the alchemical regions.
        :return:
        '''
        found_torsion = False
        forces_to_remove = []
        for force_index, force in enumerate(system.getForces()):
            if isinstance(force, mm.PeriodicTorsionForce):
                new_tor = mm.PeriodicTorsionForce()
                if found_torsion:
                    raise ValueError('Two torsions found unsure how to proceed')
                found_torsion = True
                for idx in range(force.getNumTorsions()):
                    if idx in intersect_torsions:
                        pass
                    else:
                        a, b, c, d, period, phase, k = force.getTorsionParameters(idx)
                        new_tor.addTorsion(a, b, c, d, period, phase, k)
                forces_to_remove.append(force_index)

            '''
            if isinstance(force, mm.CMAPTorsionForce):
                forces_to_remove.append(force_index)
            '''

        forces_to_remove.sort()
        forces_to_remove.reverse()
        for force_index in forces_to_remove:
            system.removeForce(force_index)
        if found_torsion:
            system.addForce(new_tor)

    def set_context_to_state(self, param_vals, context, NPT=True):
        '''
        Function to set the lambda values of the stored OpenMM contexts.

        :param param_vals: dict, containing lambda values
        :param context: the context we want to modify
        :param NPT: bool, flag to see if this is the NPT or NVT context
        :return:
        '''

        if NPT:
            for lam, val in param_vals.items():
                setattr(self.NPT_compound_state, lam, val)
            self.NPT_compound_state.apply_to_context(context)
        else:
            for lam, val in param_vals.items():
                setattr(self.NVT_compound_state, lam, val)
            self.NVT_compound_state.apply_to_context(context)


    def get_gradients(self, param, val, context, h, analitic_sterics=False):
        '''
        Function to compute the gradients of the potential w.r.t the alchemical parameters.

        :param param: string, for lambda parameter eg 'lambda_electrostatics_appear'
        :param val: float, value of lambda parameter
        :param context: OpenMM Context
        :param context: OpenMM Context
        :param h: float, finite difference to use
        :param analitic_sterics: boolean, are analytic steric gradients calculated (experimental)
        :return: float, gradient calculated with numerical finite difference
        '''

        #Sterics can use analytic derivative
        if 'sterics' in param and analitic_sterics:
            setattr(self.NPT_compound_state, param, val)
            self.NPT_compound_state.apply_to_context(context)
            dEdlam = context.getState(getParameterDerivatives=True).getEnergyParameterDerivatives()
            f_prime_x = dEdlam[param] * unit.kilojoule_per_mole

        #electrostatics must use numerical finite differnce
        else:
            #softcore is not granteed to be defined for values of lambda < 0, > 1
            if val == 1 and 'sterics' in param:
                # use backwards difference
                h_diff = h
                h = [val, val - h]
            elif val == 0 and 'sterics' in param:
                # use forward diff
                h_diff = h
                h = [val + h, val]
            else:
                # Use central diff
                h_diff = h
                h = [val + (0.5 * h), val - (0.5 * h)]

            setattr(self.NPT_compound_state, param, h[0])
            self.NPT_compound_state.apply_to_context(context)
            fx0 = context.getState(getEnergy=True).getPotentialEnergy()

            setattr(self.NPT_compound_state, param, h[1])
            self.NPT_compound_state.apply_to_context(context)
            fx1 = context.getState(getEnergy=True).getPotentialEnergy()

            # return to original param
            setattr(self.NPT_compound_state, param, val)
            self.NPT_compound_state.apply_to_context(context)

            f_prime_x = (fx0 - fx1) / h_diff
        return f_prime_x

    def build_simulation(self, system, device_id):
        '''
        Function to take OpenMM system and build OpenMM simulation object.

        :param system: OpenMM system object
        :param device_id: str, index for which cuda device this simulation will run on
        :return: dict, containing OpenMM simulation and integrator.
        '''

        integrator = mm.LangevinIntegrator(self.temp, self.friction, self.time_step)
        integrator.setConstraintTolerance(0.00001)

        platform = mm.Platform.getPlatformByName(self.platform)
        if self.platform == 'CUDA':
            properties = {'CudaPrecision': 'mixed', 'CudaDeviceIndex': device_id}
        elif self.platform == 'OpenCL':
            properties = {'OpenCLPrecision': 'mixed', 'OpenCLDeviceIndex': device_id}
        elif self.platform == 'CPU':
            properties = {}
        else:
            raise ValueError('Unknown platform {} in ties. Please select from CPU/CUDA/OpenCL'.format(self.platform))

        sim = app.Simulation(self.topology_file.topology, system, integrator, platform, properties)

        return {'sim': sim, 'integrate': integrator}

    def debug_force(self, system):
        '''
        Function which removes all but nonbonded forces while maintaining NPT ensemble

        :param system: OpenMM system to modify
        :return:
        '''
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

    def add_consraints(self):
        '''
        Function to add constraints to OpenMM system. Uses values for constrains initialized from file during class
        construction. Note: as written these constraints may not work with periodic boundary conditions on GPU.

        '''
        force = mm.CustomExternalForce("scale*k*periodicdistance(x, y, z, x0, y0, z0)^2")
        force.addGlobalParameter("scale", 0.0)
        force.addPerParticleParameter("k")
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")
        positions = self.og_positions
        force_constants = self.constraints
        if len(positions) != len(force_constants):
            raise ValueError("List of atom position and list of constraint force constants found to be different lengths")
        for i, (pos, con) in enumerate(zip(positions, force_constants)):
            force.addParticle(i, [con, *pos])
        self.NVT_alchemical_system.addForce(force)


class PDB_line(object):
    '''
    Class to let us address line of a PDB file in a readable format

    :param str: line from a pdb file

    '''
    def __init__(self, line):
        self.ah_f = str(line[0:6])
        self.atom_serial_number = int(line[6:11])
        self.atom_name = str(line[12:16])
        self.alternate_location_indicator = str(line[16:17])
        self.residue_name = str(line[17:20])
        self.chain_identifier = str(line[21:22])
        self.residue_sequence_number = int(line[22:26])
        self.code_insertion_residues = str(line[26:27])
        self.X = float(line[30:38])
        self.Y = float(line[38:46])
        self.Z = float(line[46:54])
        self.occupancy = float(line[54:60])
        try:
            self.temperature_factor = float(line[60:66])
        except ValueError:
            self.temperature_factor = None
        self.element_symbol = str(line[76:78])
        self.charge_atom = str(line[78:80])


def get_alchemical_atoms(pdb_file, absolute):
    '''
    Function to read PDB file and pull out alchemical indexes from temp factor or occupancy column (set by user).

    :param pdb_file: str, pointing to the PDB file to read
    :param absolute: bool, used in logic to determine how to read idxs

    :return: List, List of ints addressing appearing and disappearing alchemical atoms
    '''
    appear = []
    disappear = []
    with open(pdb_file) as f:
        # read file
        for line in f:
            if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':
                pdb_line_data = PDB_line(line)

                if not absolute:
                    if pdb_line_data.temperature_factor == -1.0:
                        disappear.append(pdb_line_data.atom_serial_number-1)
                    elif pdb_line_data.temperature_factor == 1.0:
                        appear.append(pdb_line_data.atom_serial_number-1)
                else:
                    if pdb_line_data.residue_name == 'MOL':
                        disappear.append(pdb_line_data.atom_serial_number-1)
    return appear, disappear


def get_constraints(file, column):
    '''
    Fuction to read constraints from PDB file.

    :param file: str, file path to pdb containing constraint information
    :param column: str, determines whether the occupancy row of temp factor will be read for data

    :return: list, list of floats for constraint strength on each atom
    '''

    constraints = []
    with open(file) as f:
        for line in f:
            if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':
                pdb_line_data = PDB_line(line)
                if column == 'occupancy':
                    constraints.append(pdb_line_data.occupancy*unit.kilocalories_per_mole/unit.angstrom**2)
                elif column == 'beta_factor':
                    if pdb_line_data.temperature_factor is None:
                        raise ValueError('No temperature factor in PDB')
                    constraints.append(pdb_line_data.temperature_factor*unit.kilocalories_per_mole/unit.angstrom**2)
                else:
                    raise ValueError('Unknown constraint columns, please chose from occupancy/beta_factor')
    return constraints

class System_ID(object):
    '''
    Class as a ID for a simulation providing information for what GPU this simulation should run on and
    what number repeat this simulation is out of some total number of repeats.

    :param device_id: int, for OpenMM GPU device id
    :param node_id: str, id number denoting which replica this is
    '''

    def __init__(self, device_id, node_id):
        self.device_id = str(device_id)
        self.node_id = str(node_id)


def minimization(NVT, constraint):
    '''
    Performs minimization and constraint relaxation.

    :param NVT: list[context, integrator] for NVT system
    :param constraint: list or None, indicates whether this system has constraints

    '''
    #Assumed that position and state is set correctly

    if constraint is not None:
        constraint_factors = [5, 4, 3, 2, 1, 0]
        for i, constraint_factor in enumerate(constraint_factors):
            NVT['sim'].context.setParameter("scale", constraint_factor)
            print('Minimization step {}: constraint scale = {}'.format(i, NVT['sim'].context.getParameter("scale")))
            mm.LocalEnergyMinimizer.minimize(NVT['sim'].context, maxIterations=1000)

        print('Minimization end: constraint scale = {}'.format(NVT['sim'].context.getParameter("scale")))
        mm.LocalEnergyMinimizer.minimize(NVT['sim'].context)

    else:
        print('Minimizing')
        mm.LocalEnergyMinimizer.minimize(NVT['sim'].context)


def equilibriation(NVT, NPT, steps, save_file, constraint):
    '''
    Perform NVT and then NPT equilibration for total time defined by user.

    :param NVT: dict, containing OpenMM simulation and integrator objects for NVT simulation
    :param NPT: dict, containing OpenMM simulation and integrator objects for NPT simulation
    :param steps: int, number of equilibration steps
    :param save_file: string, where to write the equilibrating state
    :param constraint: list or None, indicates whether this system has constraints

    '''
    # Assumed that positions of (NVT) and state and are set correctly

    # MAKE FRACTION USER DEFINED
    NVT_steps = int(steps * 0.01)
    NPT_steps = int(steps)

    ##NVT
    print('Performing NVT equilibration for {}'.format(NVT_steps*2*unit.femtosecond))
    tic = time.perf_counter()

    #initilize constraints
    if constraint is not None:
        NVT['sim'].context.setParameter("scale", 0.0)
        print('Constraint scale in NVT = {}'.format(NVT['sim'].context.getParameter("scale")))

    #NVT warming
    final_temperature = NVT['integrate'].getTemperature() / unit.kelvin
    initial_temperature = 50
    NVT['sim'].context.setVelocitiesToTemperature(initial_temperature)
    temps = np.linspace(initial_temperature, final_temperature, 10)
    print('Warming NVT from {}k to {}k'.format(initial_temperature, final_temperature))
    for temp in temps:
        NVT['integrate'].setTemperature(temp)
        NVT['integrate'].step(int(NVT_steps / len(temps)))

    tock = time.perf_counter()
    print('Took {} seconds'.format(tock - tic))

    #save positions and velocities to transfer to NPT system
    pos_vel = NVT['sim'].context.getState(getPositions=True, getVelocities=True)
    pos, vel = pos_vel.getPositions(), pos_vel.getVelocities()

    ##NPT
    #Set equilibriated pos and vel in NPT context
    NPT['sim'].context.setPositions(pos)
    NPT['sim'].context.setVelocities(vel)

    print('Performing NPT equilibration for {}'.format(NPT_steps*2*unit.femtosecond))
    tic = time.perf_counter()

    #constraints in dynamics is causing instabilities
    '''
    #relax constraints
    constraint_coeff = [1, 0.75, 0.5, 0.25, 0]
    for factor in constraint_coeff:
        NPT['sim'].context.setParameter("scale", constraint_coeff)
        print('constraint scale in NPT = {}'.format(NPT['sim'].context.getParameter("scale")))
        NPT['integrate'].step(NPT_steps/25)

    #run remaining steps
    NPT['integrate'].step(NPT_steps-(NPT_steps*25*len(constraint_coeff)))
    '''
    NPT['sim'].step(NPT_steps)

    tock = time.perf_counter()
    print('Took {} seconds'.format(tock - tic))

    print('Saving equilibrated state to disk here: {}'.format(save_file))
    NPT['sim'].saveState(save_file)


def preproduction(NVT, NPT, equili_steps, equili_state_file,  constraint):
    '''
    Meta function to wrap minimization and equilibration together.

    :param NVT: dict, containing OpenMM simulation and integrator objects for NVT simulation
    :param NPT: dict, containing OpenMM simulation and integrator objects for NPT simulation
    :param equili_steps: int, number of total NVT and NPT equilib steps to perform
    :param equili_state_file: string, where to write the equilibrating state
    :param constraint: list or None, indicates whether this system has constraints
    '''
    minimization(NVT, constraint)
    equilibriation(NVT, NPT, equili_steps, equili_state_file, constraint)

def add_simulation_reporters(sim, total_steps, save):
    '''
    Function to add reporters to a simulation.

    :param sim: OpenMM simulation object to add a reporter to
    :param total_steps: int for the number of total steps in the simulation
    :param save: string pointing to the file location to save dcd
    '''
    if save is not None:
        log_out = save+'.log'
    else:
        log_out = stdout

    sim.reporters.append(app.StateDataReporter(log_out, 5000, step=True,
                                               potentialEnergy=True, volume=True, temperature=True, progress=True,
                                               remainingTime=True, speed=True, totalSteps=total_steps, separator='\t'))

    if save is not None:
        sim.reporters.append(app.DCDReporter(save+'.dcd', 5000))

def remove_simulation_reporters(sim):
    '''

    :param sim: OpenMM simulation object to remove reporter from
    '''
    del sim.reporters
    sim.reporters = []


def simulate_system(ids, alch_sys, Lam, mask, cwd, niter, equili_steps, steps_per_iter=1000):

    '''
    Main function to perform simulations. Will iterate over windows assigned to this worker. Performs pre-production,
    and then production simulation. Will collect required grads and potentials for TI/FEP. Will then write any results
    to disk.

    Sampling per repeat = (number of windows * niter * steps_per_iter * 2fs) + (number of windows * equili_steps * 2fs)

    :param ids: System_ID class containing ids for GPU, repeat and node
    :param alch_sys: AlchSys class containing the alchemical system
    :param Lam: Lambda class containing lambda schedule
    :param mask: list containing ints for start and end range of windows to be run
    :param name: string, name of the current experiment we are running
    :param niter: int, number of iterations per window #How many times we sample grad
    :param equili_steps: int, number of 2fs steps for equilibration
    :param steps_per_iter: int, number of 2fs steps per iteration #How much sampling between sampling grad

    '''

    beta = 1.0 / (unit.BOLTZMANN_CONSTANT_kB * alch_sys.temp)

    print('Running simulation on device {}'.format(ids.device_id))
    print('Running replica {} of windows {}'.format(ids.node_id, list(range(mask[0], mask[1]))))

    nlambda = len(Lam.schedule[0])
    nstates = len(Lam.schedule[mask[0]:mask[1]])
    all_states = len(Lam.schedule)

    # prep data structures for results
    if 'TI' in alch_sys.methods:
        #Build numpy array to hold all our calculated grads
        grads = np.zeros([nstates, nlambda, niter], np.float64)
    else:
        grads = None
    if 'FEP' in alch_sys.methods:
        u_kln = np.zeros([nstates, all_states, niter], np.float64)
    else:
        u_kln = None

    # This creates two contexts on the GPU simultaneously which wastes VRAM memory but this may be better than continuously
    # creating and destroying these contexts.
    NVT = alch_sys.build_simulation(alch_sys.NVT_alchemical_system, ids.device_id)
    NPT = alch_sys.build_simulation(alch_sys.NPT_alchemical_system, ids.device_id)

    #test simulation for stability
    try:
        alch_sys.test_sim(NVT)
    except:
        #if test has faild try shifting some alchemical atoms
        alch_sys.shift_alchemical_positions()

    total_sim_NVT = int(0.01 * equili_steps) * (mask[1] - mask[0])
    total_sim_NPT = int(equili_steps) * (mask[1] - mask[0]) + steps_per_iter * niter * (mask[1] - mask[0])

    # Dynamics loops
    #Iterate over alchemical states
    for i, param_vals_i in enumerate(Lam.schedule[mask[0]:mask[1]]):

        #Init position of atoms
        NVT['sim'].context.setPositions(alch_sys.og_positions)
        NPT['sim'].context.setPositions(alch_sys.og_positions)

        #init state of context
        alch_sys.set_context_to_state(Lam.schedule[mask[0]:mask[1]][i], NVT['sim'].context, NPT=False)
        alch_sys.set_context_to_state(Lam.schedule[mask[0]:mask[1]][i], NPT['sim'].context, NPT=True)

        #init size of box, NVT never changes so dont need to init every state
        NPT['sim'].context.setPeriodicBoxVectors(*alch_sys.PBV)

        #try to load the equilibriation of this state performed on a differnt node or repeat
        #This might be dangerious as this file could be being written to?
        equili_file = 'state{}*.xml'.format(i+mask[0])
        equili_file = os.path.join(cwd, 'equilibration', equili_file)
        equili_files = list(glob.iglob(equili_file))

        load = False
        print('Load flag for equilibrated data is {}'.format(load))
        if len(equili_files) > 0 and load:
            # load the equilibriated state from disk
            print('Loading equilibriated state from disk located here: {}'.format(equili_files[0]))
            NPT['sim'].loadState(equili_files[0])
            #randomize velocities after load of equilibriate state
            #If velocities are randmozied is this still equilibriated?
            NPT['sim'].context.setVelocitiesToTemperature(alch_sys.temp)
        else:
            equili_file = os.path.join(cwd, 'LAMBDA_{}'.format(Lam.str_lams[i+mask[0]]), 'rep{}'.format(ids.node_id),
                                       'equilibration', 'state')
            equili_state_f = equili_file+'_NPT.xml'

            add_simulation_reporters(NVT['sim'], total_sim_NVT, save=equili_file+'_NVT')
            add_simulation_reporters(NPT['sim'], total_sim_NPT, save=equili_file+'_NPT')

            # minimize and equilibriate then save equilibriate state to disk
            preproduction(NVT, NPT, equili_steps, equili_state_f, alch_sys.constraints)

            remove_simulation_reporters(NVT['sim'])
            remove_simulation_reporters(NPT['sim'])

            # load the equilibriated state from disk, this sets the NPT vel and pos
            NPT['sim'].loadState(equili_state_f)

        #add reporter to simulation
        log = os.path.join(cwd, 'LAMBDA_{}'.format(Lam.str_lams[i+mask[0]]), 'rep{}'.format(ids.node_id), 'simulation',
                           'log')
        add_simulation_reporters(NPT['sim'], total_sim_NPT, save=log)

        #Production
        for iteration in range(niter):
            print('Propagating iteration {}/{} in state {}/{}'.format(iteration + 1, niter, i+mask[0], all_states))

            # propogate system in current state
            NPT['sim'].step(steps_per_iter)

            if 'TI' in alch_sys.methods:
                #Compute grads
                #get grad in OpenMM default units then convert to kcal/mol then remove explicit unit so we can do numpy math later
                grad = [alch_sys.get_gradients(param, val, NPT['sim'].context, 0.0001)
                            .in_units_of(unit.kilocalorie_per_mole)/unit.kilocalorie_per_mole
                        for param, val in param_vals_i.items()]

                grads[i, :, iteration] = grad

            if 'FEP' in alch_sys.methods:
                for j, param_vals_j in enumerate(Lam.schedule):
                    alch_sys.set_context_to_state(param_vals_j, NPT['sim'].context)
                    state = NPT['sim'].context.getState(getEnergy=True)
                    volume = state.getPeriodicBoxVolume()
                    potential = state.getPotentialEnergy() / unit.AVOGADRO_CONSTANT_NA
                    u_kln[i, j, iteration] = beta*(potential + alch_sys.pressure*volume)
                alch_sys.set_context_to_state(param_vals_i, NPT['sim'].context)

        remove_simulation_reporters(NPT['sim'])

    #Save results to disk
    for i, j in enumerate(Lam.str_lams[mask[0]: mask[1]]):
        if 'TI' in alch_sys.methods:
            file = os.path.join(cwd, 'LAMBDA_{}'.format(j), 'rep{}'.format(ids.node_id), 'results',
                                'TI.npy')
            print('Saving {} result to disk'.format(file))
            np.save(file, grads[i, :, :])

        if 'FEP' in alch_sys.methods:
            file = os.path.join(cwd, 'LAMBDA_{}'.format(j), 'rep{}'.format(ids.node_id), 'results',
                                'FEP.npy')
            print('Saving {} result to disk'.format(file))
            np.save(file, u_kln[i, :, :])

    #clean up
    del grads, u_kln
    del NPT, NVT

    return True

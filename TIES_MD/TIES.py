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

from .lambdas import Lambdas
try:
    from .alch import AlchSys, System_ID, simulate_system
    openmm_av = True
except ImportError:
    openmm_av = False

try:
    from openmm import unit, Vec3
except ImportError:  # OpenMM < 7.6
    from simtk import unit
    from simtk.openmm import Vec3

from functools import partial
from multiprocess import Pool
import numpy as np

import os
import sys
import time
from pathlib import Path

import importlib.resources as pkg_resources
#split means we will explitlly deal with reps in the submission
#for non split namd or TIES_MD will handle the parallisim of replicas
from .eng_scripts import namd_sub, namd_sub_split, openmm_sub_split

class TIES(object):
    '''
    Class to control TIES protocol, initializes variables and calls functions to start simulation or write input scripts

    :param cwd: str, for current working directory
    :param run_type: str, flag to say if we should run dynamics or not
    :param exp_name: str, for the names of experiment i.e. complex -> complex.pdb/complex.prmtop
    :param devices: list, list of ints for which cuda devices to use
    :param node_id: str, id denoting what replica of this simulation this execution of TIES_MD should run
    :param windows_mask: list containing ints for start and end range of windows to be run
    :param periodic: boolean determines if the simulation will be periodic
    :param lam: Lambda class, allow passing of custom lambda schedule
    :param platform: sting determines what platform OpenMM will target allowed values are ['CPU', 'CUDA', 'OpenCL']
    :param **kwargs: dict, containing setting from config file

    '''
    def __init__(self, cwd, exp_name, run_type='class', devices=None, node_id=None, windows_mask=None, periodic=True,
                 lam=None, platform='CUDA', **kwargs):

        nice_print('TIES')
        #check what engine we are using to test input args
        engine = kwargs['engine'].lower()

        #check all the config file args we need are present
        args_list = ['engine', 'temperature', 'pressure', 'sampling_per_window', 'equili_per_window', 'methods',
                       'total_reps', 'reps_per_exec', 'elec_edges', 'ster_edges', 'global_lambdas', 'constraint_file',
                       'constraint_column', 'input_type', 'box_type']

        optional_args = ['cell_basis_vec1', 'cell_basis_vec2', 'cell_basis_vec3', 'edge_length']

        openmm_args = []
        namd_args = ['version']

        if engine == 'namd':
            args_list.extend(namd_args)
        elif engine == 'openmm':
            args_list.extend(openmm_args)
            if not openmm_av:
                raise ImportError('OpenMMTools not installed please see:'
                                 ' https://UCL-CCS.github.io/TIES_MD/install.html#ties-openmm')

        # check we have all required arguments
        for argument in args_list:
            if argument not in kwargs.keys():
                raise ValueError('Missing option {} in configuration file'.format(argument))

        # check we have no unexpected arguments
        for argument in kwargs.keys():
            if argument not in args_list+optional_args:
                raise ValueError('Argument {} not supported for this engine or at all.'
                                 ' Please remove from the TIES.cfg.'.format(argument))

        self.all_args = args_list+optional_args

        #Iterate over our args_dict to set attributes of class to values in dict
        print('Running with arguments:')
        for k, v in kwargs.items():
            print('{} = {}'.format(k, v))
            setattr(self, k, v)

        #deal with engine common config file input
        supported_engines = ['openmm', 'namd']
        self.engine = self.engine.lower()
        if self.engine not in supported_engines:
            raise ValueError('Engine {} not supported please select from {}'.format(self.engine, supported_engines))

        self.temperature = self.temperature.split('*unit.')
        self.temperature = unit.Quantity(float(self.temperature[0]), getattr(unit, self.temperature[1]))

        self.pressure = self.pressure.split('*unit.')
        self.pressure = unit.Quantity(float(self.pressure[0]), getattr(unit, self.pressure[1]))

        self.sampling_per_window = self.sampling_per_window.split('*unit.')
        self.sampling_per_window = unit.Quantity(float(self.sampling_per_window[0]), getattr(unit, self.sampling_per_window[1]))

        self.equili_per_window = self.equili_per_window.split('*unit.')
        self.equili_per_window = unit.Quantity(float(self.equili_per_window[0]), getattr(unit, self.equili_per_window[1]))

        self.elec_edges = self.elec_edges.split(',')
        self.elec_edges = [float(x) for x in self.elec_edges]

        self.ster_edges = self.ster_edges.split(',')
        self.ster_edges = [float(x) for x in self.ster_edges]

        self.methods = self.methods.split(',')

        self.total_reps = int(self.total_reps)
        self.reps_per_exec = int(self.reps_per_exec)

        self.global_lambdas = [float(x) for x in self.global_lambdas.split(',')]
        self.windows = len(self.global_lambdas)

        if 'na' == self.constraint_file:
            self.constraint_file = None
            self.constraint_column = None

        if self.box_type == 'na':
            vecs = ['cell_basis_vec1', 'cell_basis_vec2', 'cell_basis_vec3']
            for vec in vecs:
                if vec not in kwargs.keys():
                    raise ValueError('If box type is unspecified as na in TIES.cfg the box vectors must be manually specified.'
                                     ' Please add options {} {} {} to TIES.cfg'.format(*vecs))
            self.cell_basis_vec1 = [float(x) for x in self.cell_basis_vec1.split(',')]
            self.cell_basis_vec2 = [float(x) for x in self.cell_basis_vec2.split(',')]
            self.cell_basis_vec3 = [float(x) for x in self.cell_basis_vec3.split(',')]

            self.basis_vectors = [Vec3(*self.cell_basis_vec1)*unit.angstrom,
                                  Vec3(*self.cell_basis_vec2)*unit.angstrom,
                                  Vec3(*self.cell_basis_vec3)*unit.angstrom]

        else:
            print('Getting box vectors for {} box. Ignoring cell basis vectors in TIES.cfg.'.format(self.box_type))
            if 'edge_length' not in kwargs.keys():
                raise ValueError('Must provide edge_length option in TIES.cfg to compute box vectors. If custom box vectors'
                                 ' are desired set box_type = na in TIES.cfg.')
            self.edge_length = self.edge_length.split('*unit.')
            self.edge_length = unit.Quantity(float(self.edge_length[0]), getattr(unit, self.edge_length[1]))
            self.edge_length = self.edge_length.in_units_of(unit.angstrom) / unit.angstrom
            self.basis_vectors = get_box_vectors(self.box_type, self.edge_length)

            self.cell_basis_vec1 = [float(x) for x in self.basis_vectors[0]/unit.angstrom]
            self.cell_basis_vec2 = [float(x) for x in self.basis_vectors[1]/unit.angstrom]
            self.cell_basis_vec3 = [float(x) for x in self.basis_vectors[2]/unit.angstrom]

        #deal with openmm specific config file input
        if self.engine == 'openmm':
            self.absolute = False
            #hidden for now
            #self.absolute = self.absolute.lower()
            #if self.absolute == 'yes':
            #    self.absolute = True
            #elif self.absolute == 'no':
            #    self.absolute = False
            #else:
            #    raise ValueError('Absolute option must be set to yes or no')

        else:
            self.absolute = None

        # deal with namd specific config file input
        if self.engine == 'namd':
            self.namd_version = float(self.version)
            if self.input_type != 'AMBER':
                raise ValueError('Only AMBER input supported in NAMD version of TIES MD')
        else:
            self.namd_version = None

        #Deal with command line input.
        # Input for all engines is processed here if it does not pertain to the current engine its set here but not used
        if windows_mask is None:
            self.windows_mask = [0, self.windows]
        else:
            self.windows_mask = windows_mask
        if devices is None:
            self.devices = [0]
        else:
            self.devices = devices
        self.platform = platform

        self.cwd = cwd
        self.node_id = node_id
        self.exp_name = exp_name
        self.periodic = periodic

        #build schedule for lambdas
        if lam is None:
            # build lambda schedule
            self.lam = Lambdas(self.elec_edges, self.ster_edges, self.global_lambdas)
        else:
            self.lam = lam
            print('Using custom lambda schedule all schedule options will be ignored')

        #deal with how hcp will run jobs, submission scripts etc...
        self.split_run = False
        #The user wants a total number of reps but they want to run many instances of TIES_MD each handeling 1 run
        if self.total_reps != self.reps_per_exec:
            self.split_run = True

            if self.reps_per_exec != 1:
                raise ValueError('If you wish to run a subset of repeats per execution of TIES MD please only '
                                 'use one replica per execution and set reps_per_exec=1 in TIES.cfg')

            #check the user has not given too many GPUS for the one replica
            if len(self.devices) > 1:
                raise ValueError('1 replica per execution has been specified in TIES.cfg but multiple CUDA devices'
                                 'have been specified on the command line. Please only specify 1 device.')

            #check if we are about to run OpenMM and other instances of TIES_MD could be running we have
            # node_id set such that the output is writen to a unique location
            if not self.engine == 'namd':
                if self.node_id is None and run_type == 'run':
                    raise ValueError('If total_reps != reps_per_exec then the command line option --node_id'
                                     ' must be set. Please set --node_id=X where X is an integer describing which replica '
                                     'this execution of TIES_MD should be running.')

        self.sub_header, self.sub_run_line = get_header_and_run(self.engine, self.namd_version, self.split_run,
                                                                self.global_lambdas, self.total_reps, self.exp_name)
        #Perform a job
        if run_type == 'run':
            if self.engine == 'namd':
                raise ValueError('Can only run TIES NAMD with option --run_type=setup. After initial setup jobs should'
                                 ' be submitted via the HPC scheduler. Please see example submission scripts here:'
                                 ' https://UCL-CCS.github.io/TIES_MD/HPC_submissions.html')
            TIES.run(self)
            TIES.write_analysis_cfg(self)

        elif run_type == 'setup':
            TIES.setup(self)

        elif run_type == 'class':
            print('Experiments {} initialized from dir {}'.format(self.exp_name, self.cwd))
            TIES.write_analysis_cfg(self)
        else:
            raise ValueError('Unknown run method selected from run/setup/class')

        if run_type != 'class':
            nice_print('END')

    def build_results_dirs(self, folders):
        '''
        Helper function to build output directories.
        :param folders: List of strings for what folders to build
        :return: None
        '''
        # If output folders do not exist make them.
        rep_dirs = range(self.total_reps)

        for i in rep_dirs:
            for lam in range(self.windows):
                print('Attempting to make eq, sim and results folders for LAMBDA_{}/rep{}'.format(lam, i))
                lam_dir = 'LAMBDA_{}'.format(lam)
                for folder in folders:
                    path = os.path.join(self.cwd, lam_dir, 'rep{}'.format(i), folder)
                    Path(path).mkdir(parents=True, exist_ok=True)

    def get_options(self):
        for arg in self.all_args:
            print('{}: {}'.format(arg, self.__getattribute__(arg)))

    def setup(self):
        '''
        Function to setup simulations and then stop

        '''
        if self.engine == 'namd':
            folders = ['equilibration', 'simulation']
            path = os.path.join(self.cwd, 'replica-confs')
            Path(path).mkdir(parents=True, exist_ok=True)
            self.write_namd_scripts()
        else:
            folders = ['equilibration', 'simulation', 'results']
            self.write_openmm_submission()
        TIES.build_results_dirs(self, folders)
        self.write_analysis_cfg()

    def run(self):
        '''
        Function to run the simulations in parallel

        '''
        folders = ['equilibration', 'simulation', 'results']
        TIES.build_results_dirs(self, folders)

        system = AlchSys(self.cwd, self.exp_name, self.temperature, self.pressure, self.constraint_file,
                         self.constraint_column, self.methods, self.basis_vectors, self.input_type, self.absolute,
                         self.periodic, self.platform)

        if self.split_run:
            system_ids = [System_ID(self.devices[0], self.node_id)]
        else:
            system_ids = [System_ID(self.devices[i % len(self.devices)], i) for i in range(self.total_reps)]
            
        #check output dirs exist
        for rep in system_ids:
            for lam in range(self.windows):
                lam_dir = 'LAMBDA_{}'.format(lam)
                path = os.path.join(self.cwd, lam_dir, 'rep{}'.format(rep.node_id))
                if not os.path.exists(path):
                    raise ValueError('Output dir {} missing. Command line option --node_id may be set incorrectly'.format(path))

        func = partial(simulate_system, alch_sys=system, Lam=self.lam, mask=self.windows_mask,
                       cwd=self.cwd, niter=int(self.sampling_per_window/(2.0*unit.picosecond)),
                       equili_steps=int(self.equili_per_window/(2.0*unit.femtoseconds)))
        pool = Pool(processes=len(self.devices))

        tic = time.perf_counter()
        simulated = pool.map(func, system_ids) ### MAIN MD
        toc = time.perf_counter()
        total_simulation_time = (toc - tic)

        # Kill extra threads
        pool.close()
        pool.join()
        pool.terminate()

        num_windows = self.windows_mask[1] - self.windows_mask[0]
        total_sampling = (((self.sampling_per_window +
                            self.equili_per_window) * num_windows * self.reps_per_exec) / len(self.devices)) / unit.nanosecond

        speed = total_sampling / total_simulation_time
        speed *= 86400  # seconds in a day
        print(f'Total of {total_sampling:0.1f} ns simulated in {toc - tic:0.4f} seconds')
        #The differnce between the time reported here and the time in the logs is from alchemy
        print('Simulation speed per GPU = {} ns/day'.format(speed))

        nice_print('Finished simulation')

    def write_analysis_cfg(self):
        '''
        Function to write the input configuration files for analysis.

        '''

        print('Writing analysis config files...')

        vdw_a = ','.join(str(x) for x in self.lam.lambda_sterics_appear)
        ele_a = ','.join(str(x) for x in self.lam.lambda_electrostatics_appear)
        vdw_d = ','.join(str(x) for x in self.lam.lambda_sterics_disappear)
        ele_d = ','.join(str(x) for x in self.lam.lambda_electrostatics_disappear)

        if self.engine == 'namd' and float(self.namd_version) < 3:
            eng = 'NAMD2'
        if self.engine == 'namd' and float(self.namd_version) >= 3:
            eng = 'NAMD3'
        if self.engine == 'openmm':
            eng = 'OpenMM'

        common_cfg = '''#Temperature of simulation in units of kelvin
temperature = {0}
#Directory where any output folders are writen.
output_dir = ./analysis
#Names of thermodynamic legs in simulation, corresponds to directory names in results.
legs = {1}
#Names of engines to make analysis for (NAMD2, NAMD3, OpenMM)
engines = {2}
#Directory where input data can be found, dir structure of results is fixed as standard TIES structure.
data_root = {3}
#File path pointing to experimental data, values of exp data can be blank but protein and ligand names provided
# determine what directories will be searched for results
exp_data = ./exp.dat
#Can provide comma separated list of numbers. This removes windows from analysis, zero indexed.
windows_mask = None
#str (TI,FEP) what methods perform analysis with for this engine
methods = {4}
#boolean to select if distributions of dG are calculated (0, 1)
distributions = 0
rep_convg = None
sampling_convg = None
namd_version = {5}
#comma separated list of floats for lambda schedule
vdw_a = {6}
ele_a = {7}
vdw_d = {8}
ele_d = {9}
        '''.format(self.temperature.in_units_of(unit.kelvin)/unit.kelvin, 'EDIT ME', eng, './', ','.join(self.methods),
                   self.namd_version, vdw_a, ele_a, vdw_d, ele_d)

        dummy_exp = '{\'SYSTEM NAME\': {\'LIGAND NAME\': [0.0, 0.0]}}'

        file_path = os.path.join(self.cwd, '../../../analysis.cfg')
        if not os.path.exists(file_path):
            with open(file_path, 'w') as f:
                f.write(common_cfg)

        file_path = os.path.join(self.cwd, '../../../exp.dat')
        if not os.path.exists(file_path):
            with open(file_path, 'w') as f:
                f.write(dummy_exp)

    def write_namd_scripts(self):
        '''
        Function to call all functions needed to write all namd input scripts

        '''
        self.write_namd_min()
        self.write_namd_eq()
        self.write_namd_prod()
        self.write_namd_submissions()

    def write_namd_min(self):
        '''
        Function to write eq0.conf scripts used in minimization stage

        '''
        
        if self.namd_version < 3:
            header = """#NAMD2
                     """
        else:
            header = """#NAMD3
CUDASOAintegrate off
                     """

        #check if we need constraints
        if self.constraint_file is not None:
            if self.constraint_column == 'occupancy':
                cons_col = 'O'
            else:
                cons_col = 'B'
            cons = """
constraints  on
consexp  2
consref  ../build/{}.pdb ;#need all positions
conskfile  ../build/{}
conskcol  {}
                """.format(self.exp_name, self.constraint_file, cons_col)
            run = """
set factor 10
set nall 10
set n 1
while {$n <= $nall} {
    constraintScaling $factor
    minimize 1000
    set n [expr $n + 1]
    set factor [expr $factor * 0.5]
}
constraintScaling 0
minimize 1000
                """
        else:
            cons = 'constraints  off'

            run = """
minimize 2000
            """

        #Put PBC into a dictionary to make it easier to write to file
        pbc = self.cell_basis_vec1+self.cell_basis_vec2+self.cell_basis_vec3
        vec_keys = ['cbv1', 'cbv2', 'cbv3', 'cbv4', 'cbv5', 'cbv6', 'cbv7', 'cbv8', 'cbv9']
        pbc_box = {}
        for k, v in zip(vec_keys, pbc):
            pbc_box[k] = v

        #check units on temp and make unit less for writing to file
        temp = self.temperature.in_units_of(unit.kelvin)/unit.kelvin

        #populate and write main script
        min_file = 'min.conf'
        if not self.split_run:
            min_namd_uninitialised = pkg_resources.open_text(namd_sub, min_file).read()
        else:
            min_namd_uninitialised = pkg_resources.open_text(namd_sub_split, min_file).read()
        min_namd_initialised = min_namd_uninitialised.format(structure_name=self.exp_name, constraints=cons, **pbc_box,
                                                             temp=temp, ele_start=self.elec_edges[0],
                                                             ster_end=self.ster_edges[1], header=header, run=run,
                                                             root=self.cwd)
        out_name = 'eq0.conf'
        open(os.path.join(self.cwd, './replica-confs', out_name), 'w').write(min_namd_initialised)
        if not self.split_run:
            # populate and write replica script which controls replica submissions
            min_namd_uninitialised = pkg_resources.open_text(namd_sub, 'min-replicas.conf').read()
            min_namd_initialised = min_namd_uninitialised.format(reps=self.total_reps, root=self.cwd)
            out_name = 'eq0-replicas.conf'
            open(os.path.join(self.cwd, './replica-confs', out_name), 'w').write(min_namd_initialised)

    def write_namd_eq(self):
        '''
        Function to write eq1 and eq2 namd scripts used for equilibration
        '''

        # prep common elements in eq1 and eq2
        if self.namd_version < 3:
            header = """#NAMD2
                     """
        else:
            header = """#NAMD3
CUDASOAintegrate on
                     """
        # check units on temp and pressure and make unit less for writing to file
        pressure_val = self.pressure.in_units_of(unit.bar) / unit.bar
        temp = self.temperature.in_units_of(unit.kelvin) / unit.kelvin
        steps = int(self.equili_per_window.in_units_of(unit.femtoseconds) / (2 * unit.femtoseconds))
        #round to nearest 10 so we are free to make even divisions to calculate restart frequency etc
        steps = round(steps, -1)

        #loop of eq stages 1 and 2
        for i in range(1, 3):
            #prep eq1 unique elements
            if i == 1:
                #add constraints if needed
                if self.constraint_file is not None:
                    run = """
constraintScaling 1
run {}
                        """.format(int(round(steps*0.01, -1)))
                else:
                    run = """
run {}
                        """.format(int(round(steps*0.01, -1)))
                res_freq = int(round(steps*0.01, -1)//2)
                #no barostat in eq1
                pressure = ''

            #prep eq2 unique elements
            else:
                # add constraints if needed
                if self.constraint_file is not None:
                    run = """
set factor 1
set nall 10
set n 1
while {{$n <= $nall}} {{
   constraintScaling $factor
   run {}
   set n [expr $n + 1]
   set factor [expr $factor * 0.5]
}}
constraintScaling 0
run {}
                        """.format(int(steps*0.04), int(steps*0.6))
                else:
                    run = """
run {}
                    """.format(int(steps))

                res_freq = int(steps/2)

                # Older TIES protocol uses Berendsen
                if self.namd_version < 2.14:
                    pressure = """
useGroupPressure      yes ;# needed for 2fs steps
useFlexibleCell       no  ;# no for water box, yes for membrane
useConstantArea       no  ;# no for water box, yes for membrane
BerendsenPressure                       on
BerendsenPressureTarget                 {}
BerendsenPressureCompressibility        4.57e-5
BerendsenPressureRelaxationTime         100
BerendsenPressureFreq                   2
                        """.format(pressure_val)
                else:
                    pressure = """
useGroupPressure      yes ;# needed for 2fs steps
useFlexibleCell       no  ;# no for water box, yes for membrane
useConstantArea       no  ;# no for water box, yes for membrane
langevinPiston          on             # Nose-Hoover Langevin piston pressure control
langevinPistonTarget  {}               # target pressure in bar 1atm = 1.01325bar
langevinPistonPeriod  200.0             # oscillation period in fs. correspond to pgamma T=50fs=0.05ps
langevinPistonTemp    {}                # f=1/T=20.0(pgamma)
langevinPistonDecay   100.0             # oscillation decay time. smaller value corresponds to larger random
                                       # forces and increased coupling to the Langevin temp bath.
                                       # Equal or smaller than piston period
                        """.format(pressure_val, temp)
            #add constraints if needed
            if self.constraint_file is not None:
                if self.constraint_column == 'occupancy':
                    cons_col = 'O'
                else:
                    cons_col = 'B'
                cons = """
constraints  on
consexp  2
# use the same file for the position reference and the B column
consref  ../build/{}.pdb ;#need all positions
conskfile  ../build/{}
conskcol  {}
                        """.format(self.exp_name, self.constraint_file, cons_col)
            else:
                cons = 'constraints  off'


            # read unpopulated eq file from disk
            eq_file = 'eq.conf'
            if not self.split_run:
                eq_namd_uninitialised = pkg_resources.open_text(namd_sub, eq_file).read()
            else:
                eq_namd_uninitialised = pkg_resources.open_text(namd_sub_split, eq_file).read()
            prev_output = 'eq{}'.format(i - 1)

            #populate eq file
            eq_namd_initialised = eq_namd_uninitialised.format(constraints=cons, output='eq%d' % i,
                                                               prev_output=prev_output, structure_name=self.exp_name,
                                                               pressure=pressure, run=run, temp=temp,
                                                               ele_start=self.elec_edges[0], ster_end=self.ster_edges[1],
                                                               header=header, res_freq=res_freq, root=self.cwd)
            out_name = "eq{}.conf".format(i)
            open(os.path.join(self.cwd, './replica-confs', out_name), 'w').write(eq_namd_initialised)

            if not self.split_run:
                #read and write eq replica to handle replica simulations
                eq_namd_uninitialised = pkg_resources.open_text(namd_sub, 'eq-replicas.conf').read()
                eq_namd_initialised = eq_namd_uninitialised.format(reps=self.total_reps,
                                                                   prev=prev_output, current='eq{}'.format(i),
                                                                   root=self.cwd)
                open(os.path.join(self.cwd, './replica-confs', 'eq{}-replicas.conf'.format(i)), 'w').write(eq_namd_initialised)

    def write_namd_prod(self):
        '''
        Function to write sim1 namd conf script for production simulation
        '''

        # check units on steps, temp and pressure and make unit less for writing to file
        pressure_val = self.pressure.in_units_of(unit.bar) / unit.bar
        temp = self.temperature.in_units_of(unit.kelvin) / unit.kelvin
        steps = int(self.sampling_per_window.in_units_of(unit.femtoseconds)/(2*unit.femtoseconds))
        # round to nearest 10
        steps = round(steps, -1)

        #set header for NAMD2 or 3
        if self.namd_version < 3:
            header = """#NAMD2
                                        """
        else:
            header = """#NAMD3
CUDASOAintegrate on                                 """

        #set barostat for NAMD2/NAMD3
        # NAMD3 cant use Berendsen
        if self.namd_version < 3:
            pressure = """BerendsenPressure                       on
BerendsenPressureTarget                 {}
BerendsenPressureCompressibility        4.57e-5
BerendsenPressureRelaxationTime         100
BerendsenPressureFreq                   2
                            """.format(pressure_val)
        else:
            pressure = """langevinPiston          on             # Nose-Hoover Langevin piston pressure control
langevinPistonTarget  {}               # target pressure in bar 1atm = 1.01325bar
langevinPistonPeriod  200.0            # oscillation period in fs. correspond to pgamma T=50fs=0.05ps
langevinPistonTemp    {}              # f=1/T=20.0(pgamma)
langevinPistonDecay   100.0            # oscillation decay time. smaller value corresponds to larger random
                                       # forces and increased coupling to the Langevin temp bath.
                                       # Equal or smaller than piston period
                            """.format(pressure_val, temp)

        # read unpopulated eq file from disk
        sim_file = 'sim1.conf'
        if not self.split_run:
            sim_namd_uninitialised = pkg_resources.open_text(namd_sub, sim_file).read()
        else:
            sim_namd_uninitialised = pkg_resources.open_text(namd_sub_split, sim_file).read()
        sim_namd_initialised = sim_namd_uninitialised.format(structure_name=self.exp_name, temp=temp, pressure=pressure,
                                                             ele_start=self.elec_edges[0], ster_end=self.ster_edges[1],
                                                             header=header, steps=steps, root=self.cwd)
        out_name = "sim1.conf"
        open(os.path.join(self.cwd, './replica-confs', out_name), 'w').write(sim_namd_initialised)

        # read and write sim replica to handle replica simulations, only if we want to use +replicas option
        if not self.split_run:
            sim_namd_uninitialised = pkg_resources.open_text(namd_sub, 'sim1-replicas.conf').read()
            sim_namd_initialised = sim_namd_uninitialised.format(reps=self.total_reps, root=self.cwd)
            open(os.path.join(self.cwd, './replica-confs', 'sim1-replicas.conf'), 'w').write(sim_namd_initialised)

    def write_namd_submissions(self):
        '''
        Function to write an example submission script of NAMD job on HPC (ARCHER2)
        '''
        lambs = [str(x) for x in self.global_lambdas]
        lambs = ' '.join(lambs)

        if self.namd_version < 3:
            if not self.split_run:
                namd_uninitialised = pkg_resources.open_text(namd_sub, 'sub.sh').read()
                namd_initialised = namd_uninitialised.format(lambs=lambs, reps=self.total_reps, run_line=self.sub_run_line,
                                                             header=self.sub_header, root=self.cwd)
                open(os.path.join(self.cwd, 'sub.sh'), 'w').write(namd_initialised)
            else:
                namd_uninitialised = pkg_resources.open_text(namd_sub_split, 'sub.sh').read()
                namd_initialised = namd_uninitialised.format(lambs=lambs, reps=self.total_reps, run_line=self.sub_run_line,
                                                             header=self.sub_header, root=self.cwd)
                open(os.path.join(self.cwd, 'sub.sh'), 'w').write(namd_initialised)

        else:
            print('No NAMD3 example script currently implemented. Examples for GPU scripts can be found here '
                  'https://UCL-CCS.github.io/TIES_MD/HPC_submissions.html')

    def write_openmm_submission(self):
        '''
        Function to write an example submission script of OpenMM job on HPC (Summit)

        :return:
        '''

        lambs = [str(x) for x in self.global_lambdas]
        lambs = ' '.join(lambs)

        if self.split_run:
            openmm_uninitialised = pkg_resources.open_text(openmm_sub_split, 'sub.sh').read()
            openmm_initialised = openmm_uninitialised.format(lambs=lambs, reps=self.total_reps, run_line=self.sub_run_line,
                                                         header=self.sub_header, root=self.cwd, py_bin=sys.path[0])
            open(os.path.join(self.cwd, 'sub.sh'), 'w').write(openmm_initialised)
        else:
            pass


def get_box_vectors(box_type, d):
    '''
    Function to compute the box vectors if we know the box type and edge length.

    :param box_type: str, defining the bix type to make vectors for
    :param d: float for edge length of box unit angstrom

    :return: list of Vec3 for each basis vector

    '''

    supported = ['cube', 'truncatedOctahedron', 'rhombicDodecahedron']
    if box_type == 'cube':
        vectors = Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)
    elif box_type == 'truncatedOctahedron':
        vectors = Vec3(1, 0, 0), Vec3(1 / 3, 2 * np.sqrt(2) / 3, 0), Vec3(-1 / 3, np.sqrt(2) / 3, np.sqrt(6) / 3)
    elif box_type == 'rhombicDodecahedron':
        vectors = Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0.5, 0.5, np.sqrt(2) / 2)
    else:
        raise ValueError('Unknown box type, please select from {}'
                         ' or specify vectors manually in TIES.cfg'.format(supported))

    return [d * v * unit.angstrom for v in vectors]


def sort_results(string):
    '''
    Helper function to sort directory paths by specific index in file name.

    :param string: Path to a results file.

    '''
    return int(string.split('/')[-4].split('_')[1])


def nice_print(string):
    '''
    Function to print str with padding.

    :param string: A string which should be printed.

    '''
    string = string.center(75, '#')
    print(string)


def get_header_and_run(engine, namd_version, split_run, global_lambdas, reps, exp_name):
    '''
    Function to prep submission file. Number of windows and replicas are inspected to make best guess at
    number of nodes and CPUS/GPUS

    :param engine: str, What engine are we using [namd, openmm]
    :param namd_version: float, What version of namd are we using if any
    :param split_run: bool, Should each run line run all or one replica
    :param global_lambdas: list of floats for value of global lambda in eah window
    :param reps: int, number of replica simulations
    :param exp_name str, name of the experiment e.g. complex
    :return: str, str for header and run line of HPC sub script
    '''

    num_windows = len(global_lambdas)

    if engine == 'namd':
        if namd_version < 3:
            #ARCHER2 specific
            num_cpu = 128
            if split_run:
                sub_header = """#Example script for ARCHER2 NAMD2
#SBATCH --job-name=LIGPAIR
#SBATCH --nodes={0}
#SBATCH --tasks-per-node={1}
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00:00
#SBATCH --account=XXX
#SBATCH --partition=standard
#SBATCH --qos=standard

module load namd/2.14-nosmp

#--nodes and nodes_per_namd can be scaled up for large simulations
nodes_per_namd=1
cpus_per_namd={1}""".format(int(num_windows*reps), num_cpu)

                sub_run_line = 'srun -N $nodes_per_namd -n $cpus_per_namd --distribution=block:block' \
                                ' --hint=nomultithread namd2 --tclmain eq$stage.conf $lambda $win_id $i&'

            else:
                sub_header = """#Example script for ARCHER2 NAMD2
#SBATCH --job-name=LIGPAIR
#SBATCH --nodes={}
#SBATCH --tasks-per-node={}
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00:00
#SBATCH --account=XXX
#SBATCH --partition=standard
#SBATCH --qos=standard

module load namd/2.14-nosmp

#--nodes and nodes_per_namd can be scaled up for large simulations
nodes_per_namd={}
cpus_per_namd={}""".format(int(num_windows*reps), num_cpu, reps, int(reps*num_cpu))

                sub_run_line = 'srun -N $nodes_per_namd -n $cpus_per_namd --distribution=block:block' \
                                ' --hint=nomultithread namd2 +replicas {} --tclmain eq$stage-replicas.conf $lambda $win_id&'.format(reps)

        else:
            #no NAMD3
            sub_header = None
            sub_run_line = None

    else:
        if split_run:
            #summit specific
            gpus_per_node = 6
            num_jobs = reps*num_windows

            sub_header = """#Example script for Summit OpenMM
#BSUB -P XXX
#BSUB -W 240
#BSUB -nnodes {}
#BSUB -alloc_flags "gpudefault smt1"
#BSUB -J LIGPAIR
#BSUB -o oLIGPAIR.%J
#BSUB -e eLIGPAIR.%J""".format(int(np.ceil(num_jobs/gpus_per_node)))
            sub_run_line = 'jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 TIES_MD --config_file=$ties_dir/TIES.cfg' \
                           ' --exp_name={} --windows_mask=$lambda,$(expr $lambda + 1)' \
                           ' --node_id=$i > $ties_dir/$lambda_$i.out&'.format(exp_name)
        else:
            # no OpenMM unified job on HPC
            sub_header = None
            sub_run_line = None

    return sub_header, sub_run_line



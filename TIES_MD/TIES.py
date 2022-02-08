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

from functools import partial
from simtk import unit
from simtk.openmm import Vec3
from multiprocess import Pool
import numpy as np

import os
import time
from pathlib import Path

import importlib.resources as pkg_resources
from .eng_scripts import namd

class TIES(object):
    '''
    Class to control TIES protocol, initializes variables and calls functions to start simulation or write input scripts

    :param cwd: str, for current working directory
    :param run_type: str, flag to say if we should run dynamics or not
    :param exp_name: str, for the names of experiment i.e. complex -> complex.pdb/complex.prmtop
    :param devices: list, list of ints for which cuda devices to use
    :param node_id: str, id denoting what replica of this simulation this execution of TIES_MD should run
    :param windows_mask: list of ints for ids for what windows to run
    :param periodic: boolean determines if the simulation will be periodic
    :param args_dict: dict, containing setting from config file
    :param lam: Lambda class, allow passing of custom lambda schedule

    '''
    def __init__(self, cwd, run_type, exp_name, devices, node_id, windows_mask, periodic, args_dict, lam=None):

        nice_print('TIES')

        #check what engine we are using to test input args
        engine = args_dict['engine'].lower()

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
            if argument not in args_dict:
                raise ValueError('Missing option {} in configuration file'.format(argument))

        # check we have no unexpected arguments
        for argument in args_dict.keys():
            if argument not in args_list+optional_args:
                raise ValueError('Argument {} not supported for this engine or at all.'
                                 ' Please remove from the TIES.cfg.'.format(argument))

        #Iterate over our args_dict to set attributes of class to values in dict
        print('Running with arguments:')
        for k, v in args_dict.items():
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
                if vec not in args_dict.keys():
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
            if 'edge_length' not in args_dict.keys():
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
        else:
            self.namd_version = None

        #Deal with command line input.
        # Input for all engines is processed here if it does not pertain to the current engine its set here but not used
        if windows_mask is None:
            self.windows_mask = [0, self.windows]
        else:
            self.windows_mask = windows_mask

        self.cwd = cwd
        self.devices = devices
        self.node_id = node_id
        self.exp_name = exp_name
        self.periodic = periodic

        self.split_run = False
        #The user wants a total number of reps but they want to run many instances of TIES_MD each handeling 1 run
        if self.total_reps != self.reps_per_exec:
            self.split_run = True
            if self.engine == 'namd':
                print('NAMD only supports total_reps = reps_per_exec please set these values equal in TIES.cfg')
                print('Will proceed with reps_per_exec = {}'.format(self.total_reps))
                self.reps_per_exec = self.total_reps
                self.split_run = False
            else:
                if self.node_id is None:
                    raise ValueError('If total_reps != reps_per_exec then the command line option --node_id'
                                     ' must be set. Please set --node_id=X where X is an integer describing which replica '
                                     'this execution of TIES_MD should be running.')
                if self.reps_per_exec != 1:
                    raise ValueError('If you wish to run a subset of repeats per execution of TIES MD please only '
                                     'use one replica per execution and set reps_per_exec=1 in TIES.cfg')
                if len(self.devices) > 1:
                    raise ValueError('1 replica per execution has been specified in TIES.cfg but multiple CUDA devices'
                                     'have been specified on the command line. Please only specify 1 device.')
        if len(self.devices) != self.reps_per_exec:
            print('For best parallel performance you may wish to consider allocating the same number of GPUS as'
                  ' requested replicas. GPUS allocated = {}, replicas requested {}'.format(len(self.devices), self.reps_per_exec))

        #build schedual for lambdas
        if lam is None:
            # build lambda schedule
            self.lam = Lambdas(self.elec_edges, self.ster_edges, self.global_lambdas)
        else:
            self.lam = lam
            print('Using custom lambda schedule all schedule options will be ignored')

        #Perform a job
        if run_type == 'run':
            if self.engine == 'namd':
                raise ValueError('Can only run TIES NAMD with option --run_type=setup. After initial setup jobs should'
                                 ' be submitted via the HPC scheduler. Please see example submission scripts here:'
                                 ' https://UCL-CCS.github.io/TIES_MD/HPC_submissions.html')
            TIES.run(self)
           
        elif run_type == 'setup':
            if self.engine == 'namd':
                folders = ['equilibration', 'simulation']
                path = os.path.join(self.cwd, 'replica-confs')
                Path(path).mkdir(parents=True, exist_ok=True)
                self.write_namd_scripts()
            else:
                folders = ['equilibration', 'simulation', 'results']
            TIES.build_results_dirs(self, folders)
        elif run_type == 'class':
            print('Experiments {} initialized from dir {}'.format(self.exp_name, self.cwd))
        else:
            raise ValueError('Unknown run method selected from run/class')
        
        if run_type != 'class':
            self.write_analysis_cfg()

        nice_print('END')

    def build_results_dirs(self, folders):
        '''
        Helper function to build output directories.
        :param folders: List of strings for what folders to build
        :return: None
        '''
        # If output folders do not exist make them.
        if self.split_run:
            rep_dirs = [self.node_id]
        else:
            rep_dirs = range(self.total_reps)

        for i in rep_dirs:
            for lam in range(self.windows):
                print('Attempting to make eq, sim and results folders for LAMBDA_{}/rep{}'.format(lam, i))
                lam_dir = 'LAMBDA_{}'.format(lam)
                for folder in folders:
                    path = os.path.join(self.cwd, lam_dir, 'rep{}'.format(i), folder)
                    Path(path).mkdir(parents=True, exist_ok=True)

    def run(self):
        '''
        Function to run the simulations in parallel

        '''
        folders = ['equilibration', 'simulation', 'results']
        TIES.build_results_dirs(self, folders)

        system = AlchSys(self.cwd, self.exp_name, self.temperature, self.pressure, self.constraint_file,
                         self.constraint_column, self.methods, self.basis_vectors, self.input_type, self.absolute,
                         self.periodic)

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

        namd_only = '''#float for namd version used
namd_version = {0}
#comma separated list of floats for lambda schedule
vdw_a = {1}
ele_a = {2}
vdw_d = {3}
ele_d = {4}
        '''.format(self.namd_version, vdw_a, ele_a, vdw_d, ele_d)

        openmm_only = '''#(True or False) if true all replicas are combined into one long time series
fep_combine_reps = False
#comma separated list of floats for lambda schedule
vdw_a = {0}
ele_a = {1}
vdw_d = {2}
ele_d = {3}
        '''.format(vdw_a, ele_a, vdw_d, ele_d)


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
#boolean to select if distributions of dG are calculated (not implemented)
distributions = False
        '''.format(self.temperature.in_units_of(unit.kelvin)/unit.kelvin, 'EDIT ME', eng, './', ','.join(self.methods))

        dummy_exp = '{\'SYSTEM NAME\': {\'LIGAND NAME\': [0.0, 0.0]}}'

        #write files
        if self.engine == 'namd':
            file_path = os.path.join(self.cwd, '../../../namd.cfg')
            if not os.path.exists(file_path):
                with open(file_path, 'w') as f:
                    f.write(namd_only)

        if self.engine == 'openmm':
            file_path = os.path.join(self.cwd, '../../../openmm.cfg')
            if not os.path.exists(file_path):
                with open(file_path, 'w') as f:
                    f.write(openmm_only)

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
        min_namd_uninitialised = pkg_resources.open_text(namd, min_file).read()
        min_namd_initialised = min_namd_uninitialised.format(structure_name=self.exp_name, constraints=cons, **pbc_box,
                                                             temp=temp, ele_start=self.elec_edges[0],
                                                             ster_end=self.ster_edges[1], header=header, run=run)
        out_name = 'eq0.conf'
        open(os.path.join('./replica-confs', out_name), 'w').write(min_namd_initialised)

        #populate and write replica script which controls replica submissions
        min_namd_uninitialised = pkg_resources.open_text(namd, 'min-replicas.conf').read()
        min_namd_initialised = min_namd_uninitialised.format(reps=self.total_reps)
        out_name = 'eq0-replicas.conf'
        open(os.path.join('./replica-confs', out_name), 'w').write(min_namd_initialised)

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

                # NAMD3 cant use Berendsen
                if self.namd_version <= 3:
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
langevinPistonPeriod  50.0             # oscillation period in fs. correspond to pgamma T=50fs=0.05ps
langevinPistonTemp    300              # f=1/T=20.0(pgamma)
langevinPistonDecay   25.0             # oscillation decay time. smaller value corresponds to larger random
                                       # forces and increased coupling to the Langevin temp bath.
                                       # Equal or smaller than piston period
                        """.format(pressure_val)
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

            # read unpopulated eq file from dis
            eq_file = 'eq.conf'
            eq_namd_uninitialised = pkg_resources.open_text(namd, eq_file).read()

            prev_output = 'eq{}'.format(i - 1)

            #populate eq file
            eq_namd_initialised = eq_namd_uninitialised.format(constraints=cons, output='eq%d' % i,
                                                               prev_output=prev_output, structure_name=self.exp_name,
                                                               pressure=pressure, run=run, temp=temp,
                                                               ele_start=self.elec_edges[0], ster_end=self.ster_edges[1],
                                                               header=header, res_freq=res_freq)

            out_name = "eq{}.conf".format(i)
            open(os.path.join('./replica-confs', out_name), 'w').write(eq_namd_initialised)


            #read and write eq replica to handle replica simulations
            eq_namd_uninitialised = pkg_resources.open_text(namd, 'eq-replicas.conf').read()
            eq_namd_initialised = eq_namd_uninitialised.format(reps=self.total_reps,
                                                               prev=prev_output, current='eq{}'.format(i))
            open(os.path.join('./replica-confs', 'eq{}-replicas.conf'.format(i)), 'w').write(eq_namd_initialised)

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

        # read unpopulated eq file from dis
        sim_file = 'sim1.conf'
        sim_namd_uninitialised = pkg_resources.open_text(namd, sim_file).read()

        #set jeader for NAMD2 or 3
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
langevinPistonPeriod  50.0             # oscillation period in fs. correspond to pgamma T=50fs=0.05ps
langevinPistonTemp    300              # f=1/T=20.0(pgamma)
langevinPistonDecay   25.0             # oscillation decay time. smaller value corresponds to larger random
                                       # forces and increased coupling to the Langevin temp bath.
                                       # Equal or smaller than piston period
                            """.format(pressure_val)

        sim_namd_initialised = sim_namd_uninitialised.format(structure_name=self.exp_name, temp=temp, pressure=pressure,
                                                             ele_start=self.elec_edges[0], ster_end=self.ster_edges[1],
                                                             header=header, steps=steps)
        out_name = "sim1.conf"
        open(os.path.join('./replica-confs', out_name), 'w').write(sim_namd_initialised)

        # read and write eq replica to handle replica simulations
        # read unpopulated eq file from dis
        sim_namd_uninitialised = pkg_resources.open_text(namd, 'sim1-replicas.conf').read()
        sim_namd_initialised = sim_namd_uninitialised.format(reps=self.total_reps)
        open(os.path.join('./replica-confs', 'sim1-replicas.conf'), 'w').write(sim_namd_initialised)

    def write_namd_submissions(self):
        '''
        Function to write an example submission script of NAMD job on HPC (SuperMUC-NG)
        '''
        lambs = [str(x) for x in self.global_lambdas]
        lambs = ' '.join(lambs)

        if self.namd_version < 3:
            namd_uninitialised = pkg_resources.open_text(namd, 'sub.sh').read()
            namd_initialised = namd_uninitialised.format(lambs=lambs, reps=self.total_reps, nodes=len(self.global_lambdas))
            open(os.path.join('./', 'sub.sh'), 'w').write(namd_initialised)

        else:
            print('No NAMD3 example script currently implemented. Examples for GPU scripts can be found here '
                  'https://UCL-CCS.github.io/TIES_MD/HPC_submissions.html')


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



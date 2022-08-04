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
from .alch import AlchSys, System_ID, simulate_system

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
import glob
import shutil

import importlib.resources as pkg_resources
#split means we will explicitly deal with reps in the submission
#for non split namd or TIES_MD will handle the parallisim of replicas
from .eng_scripts import namd_sub, namd_sub_split, openmm_sub_split, cfg_scripts

class TIES(object):
    '''
    Class to control TIES protocol, initializes variables and calls functions to start simulation or write input scripts

    :param cwd: str, for current working directory
    :param run_type: str, flag to say if we should run dynamics or not
    :param exp_name: str, for the names of experiment i.e. complex -> complex.pdb/complex.prmtop
    :param devices: list, list of ints for which cuda devices to use
    :param node_id: float, id denoting what replica of this simulation this execution of TIES_MD should run
    :param windows_mask: list containing ints for start and end range of windows to be run
    :param periodic: boolean determines if the simulation will be periodic
    :param lam: Lambda class, allow passing of custom lambda schedule
    :param platform: sting determines what platform OpenMM will target allowed values are ['CPU', 'CUDA', 'OpenCL']
    :param **kwargs: dict, containing setting from config file

    '''
    def __init__(self, cwd, exp_name, run_type='class', devices=None, node_id=None, windows_mask=None, periodic=True,
                 lam=None, platform='CUDA', **kwargs):
        nice_print('TIES')
        print('Wade, A.D., et al. 2022. Alchemical Free Energy Estimators and Molecular Dynamics Engines:'
              ' Accuracy, Precision, and Reproducibility. Journal of chemical theory and computation, 18(6), pp.3972-3987.\n')

        print('If you use OpenMM with this software please cite:')
        print('Chodera, J., et al. 2022. choderalab/openmmtools: 0.21.4. Zenodo. doi: 10.5281/zenodo.6625266.')
        print('Eastman, P., et al. 2017. OpenMM 7: Rapid development of high performance algorithms'
              ' for molecular dynamics. PLoS computational biology, 13(7), p.e1005659.\n')

        #check all the config file args we need are present
        args_list = ['engine', 'temperature', 'pressure', 'sampling_per_window', 'equili_per_window', 'methods',
                       'total_reps', 'reps_per_exec', 'elec_edges', 'ster_edges', 'global_lambdas', 'constraint_file',
                       'constraint_column', 'input_type', 'box_type']

        optional_args = ['cell_basis_vec1', 'cell_basis_vec2', 'cell_basis_vec3', 'edge_length']

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

        #engine must be delt with first to set namd_version which other options may need.
        api_sensitive = ['engine', 'total_reps', 'reps_per_exec', 'elec_edges', 'ster_edges', 'global_lambdas',
                         'box_type']

        #Iterate over our args_dict to set attributes of class to values in dict
        print('Read arguments from file:')
        for k, v in kwargs.items():
            if k in api_sensitive:
                full_k = '_'+k
            else:
                full_k = k
            print('{} = {}'.format(k, v))
            setattr(self, full_k, v)
        print('')

        #set any nonexistant optional args to None
        for option in optional_args:
            if option not in kwargs.keys():
                setattr(self, option, None)

        #set any attr the api needs
        if self._total_reps != self._reps_per_exec:
            self._split_run = True
        else:
            self._split_run = False
        self._exp_name = exp_name

        self._elec_edges = self._elec_edges.split(',')
        self._elec_edges = [float(x) for x in self._elec_edges]

        self._ster_edges = self._ster_edges.split(',')
        self._ster_edges = [float(x) for x in self._ster_edges]

        self._total_reps = int(self._total_reps)
        self._reps_per_exec = int(self._reps_per_exec)

        self._global_lambdas = [round(float(x), 2) for x in self._global_lambdas.split(',')]

        #set genral args
        self.windows = len(self.global_lambdas)
        self.temperature = self.temperature.split('*unit.')
        self.temperature = unit.Quantity(float(self.temperature[0]), getattr(unit, self.temperature[1]))

        self.pressure = self.pressure.split('*unit.')
        self.pressure = unit.Quantity(float(self.pressure[0]), getattr(unit, self.pressure[1]))

        self.sampling_per_window = self.sampling_per_window.split('*unit.')
        self.sampling_per_window = unit.Quantity(float(self.sampling_per_window[0]), getattr(unit, self.sampling_per_window[1]))

        self.equili_per_window = self.equili_per_window.split('*unit.')
        self.equili_per_window = unit.Quantity(float(self.equili_per_window[0]), getattr(unit, self.equili_per_window[1]))

        self.run_type = run_type
        self.methods = self.methods.split(',')

        if 'na' == self.constraint_file:
            self.constraint_file = None
            self.constraint_column = None

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
        self.periodic = periodic

        #run through api logic
        for prop in api_sensitive:
            setattr(self, prop, self.__getattribute__('_'+prop))

        #make sub file now all args populated
        self.sub_header, self.sub_run_line = get_header_and_run(self.engine, self.namd_version, self.split_run,
                                                                self.global_lambdas, self.total_reps, self.exp_name)

        # build schedule for lambdas do this last so passed lam can overwrite if desired
        print('Lambda schedule:')
        if lam is None:
            # build lambda schedule
            self.lam = Lambdas(self.elec_edges, self._ster_edges, self.global_lambdas)
        else:
            self.lam = lam
            print('Using custom lambda schedule all schedule options will be ignored')
        print('')

        #Perform a job
        if self.run_type == 'run':
            if self.engine == 'namd':
                raise ValueError('Can only run TIES NAMD with option --run_type=setup. After initial setup jobs should'
                                 ' be submitted via the HPC scheduler. Please see example submission scripts here:'
                                 ' https://UCL-CCS.github.io/TIES_MD/HPC_submissions.html')
            TIES.run(self)
            TIES.write_analysis_cfg(self)

        elif self.run_type == 'setup':
            TIES.setup(self)

        elif self.run_type == 'class':
            print('Experiments {} initialized from dir {}'.format(self.exp_name, self.cwd))
            TIES.write_analysis_cfg(self)
        else:
            raise ValueError('Unknown run method selected from run/setup/class')

        if self.run_type != 'class':
            nice_print('END')

    @property
    def box_type(self):
        """
        What type of simulation box is being used (cube, truncatedOctahedron, rhombicDodecahedron or na for manual)

        :return: string for box type.
        """
        return self._box_type

    @property
    def exp_name(self):
        '''
        What is the name of the experiment, this is the prefix expected on input files e.g. experiment_name.pdb

        :return: string for the experiment name.
        '''
        return self._exp_name

    @property
    def engine(self):
        '''
        What molecular dynamics engine is going to be used (openmm/namd2.14/namd3)

        :return: string for molecular dynamics engine.
        '''
        return self._engine

    @property
    def global_lambdas(self):
        '''
        The global values of lambda in the lambda schedule.

        :return: list of floats for lambda schedule
        '''
        return self._global_lambdas

    @property
    def elec_edges(self):
        '''
        here in lambda schedule (0->1) should the electrostatic potentials begin, stop appearing.

        :return: list of floats for when the electrostatic potentials begin, stop appearing.
        '''
        return self._elec_edges

    @property
    def ster_edges(self):
        '''
        Where in lambda schedule (0->1) should the Lennard_Jones potentials begin, stop appearing.

        :return: list of floats for when the Lennard_Jones potentials begin, stop appearing.
        '''
        return self._ster_edges

    @property
    def split_run(self):
        '''
        boolean that sets if each execution of TIESMD or NAMD runs all replicas or a subset.

        :return: bool for if run is split
        '''
        return self._split_run

    @property
    def total_reps(self):
        '''
        What is the total number of replicas we expect to run

        :return: int for total number of replicas.
        '''
        return self._total_reps

    @property
    def reps_per_exec(self):
        '''
        How many replicas will each run of the program execute.

        :return: int for replicas per run.
        '''
        return self._reps_per_exec

    @box_type.setter
    def box_type(self, value):
        '''
        Setting function for box type, will build manual box from cell basis vectors if user passes box type na
        :param value: str, for what box type we want

        :return: None
        '''
        self._box_type = value
        if self._box_type == 'na':
            vecs = ['cell_basis_vec1', 'cell_basis_vec2', 'cell_basis_vec3']
            for vec in vecs:
                if self.__getattribute__(vec) is None:
                    raise ValueError(
                        'If box type is unspecified as na in TIES.cfg the box vectors must be manually specified.'
                        ' Please add options {} {} {} to TIES.cfg'.format(*vecs))
            self.cell_basis_vec1 = [float(x) for x in self.cell_basis_vec1.split(',')]
            self.cell_basis_vec2 = [float(x) for x in self.cell_basis_vec2.split(',')]
            self.cell_basis_vec3 = [float(x) for x in self.cell_basis_vec3.split(',')]

            self.basis_vectors = [Vec3(*self.cell_basis_vec1) * unit.angstrom,
                                  Vec3(*self.cell_basis_vec2) * unit.angstrom,
                                  Vec3(*self.cell_basis_vec3) * unit.angstrom]

        else:
            print('Getting box vectors for {} box. Ignoring cell basis vectors'.format(self.box_type))
            if self.edge_length is None:
                raise ValueError('Must provide edge_length option in TIES.cfg to compute box vectors. If custom box vectors'
                                 ' are desired set box_type = na in TIES.cfg.')
            self.edge_length = self.edge_length.split('*unit.')
            self.edge_length = unit.Quantity(float(self.edge_length[0]), getattr(unit, self.edge_length[1]))
            self.edge_length = self.edge_length.in_units_of(unit.angstrom) / unit.angstrom
            self.basis_vectors = get_box_vectors(self.box_type, self.edge_length)

            self.cell_basis_vec1 = [float(x) for x in self.basis_vectors[0] / unit.angstrom]
            self.cell_basis_vec2 = [float(x) for x in self.basis_vectors[1] / unit.angstrom]
            self.cell_basis_vec3 = [float(x) for x in self.basis_vectors[2] / unit.angstrom]

    @exp_name.setter
    def exp_name(self, value):
        '''
        Setter for exp_name, will rebuild the submission file to reflect updated name.
        :param value: str, for the exp_name

        :return: None
        '''
        self._exp_name = value
        self.sub_header, self.sub_run_line = get_header_and_run(self._engine, self.namd_version, self._split_run,
                                                                self._global_lambdas, self._total_reps, self._exp_name)

    @engine.setter
    def engine(self, value):
        '''
        Setter for the engine type, check logic on selection and updates submission scripts
        :param value: str, for what engine we want.

        :return: None
        '''
        self._engine = value.lower()
        if self._engine == 'openmm':
            self.namd_version = None
            self.absolute = False
        else:
            self.absolute = None

        if self._engine == 'namd':
            raise ValueError('Must pass version with engine i.e namd2.14 or namd3')

        if 'namd' in self._engine:
            self.namd_version = float(self._engine.split('namd')[1])
            self._engine = 'namd'
            if self.namd_version != 2.14 or self.namd_version != 3:
                raise ValueError('Supported NAMD versions are: namd2.14/namd3')
            if self.input_type != 'AMBER':
                raise ValueError('Only AMBER input supported in NAMD version of TIES MD')

        supported_engines = ['openmm', 'namd']
        if self._engine not in supported_engines:
            raise ValueError('Engine {} not supported please select from {}'.format(self._engine, supported_engines))

        self.sub_header, self.sub_run_line = get_header_and_run(self._engine, self.namd_version, self._split_run,
                                                                self._global_lambdas, self._total_reps, self._exp_name)

    @ster_edges.setter
    def ster_edges(self, value):
        '''
        Setter for ster_edges, rebuilds Lambdas class with updated schedule.
        :param value: list of floats, for where the Lennard_Jones potentials begin, stop appearing.

        :return: None
        '''
        self._ster_edges = value
        self.lam = Lambdas(self._elec_edges, self._ster_edges, self._global_lambdas, debug=False)

    @elec_edges.setter
    def elec_edges(self, value):
        '''
        Setter for elec_edges, rebuilds Lambdas class with updated schedule.
        :param value: list of floats, for where the electrostatic potentials begin, stop appearing.

        :return: None
        '''
        self._elec_edges = value
        self.lam = Lambdas(self._elec_edges, self._ster_edges, self._global_lambdas, debug=False)

    @global_lambdas.setter
    def global_lambdas(self, value):
        '''
        Setter for global_lambdas, rebuilds Lambdas class and submission scripts with updated info.
        :param value: list of floats for the value the global controlling parameter takes in each window.

        :return: None
        '''
        self._global_lambdas = [round(x, 2) for x in value]
        self.lam = Lambdas(self._elec_edges, self._ster_edges, self._global_lambdas, debug=False)
        self.sub_header, self.sub_run_line = get_header_and_run(self._engine, self.namd_version, self._split_run,
                                                                self._global_lambdas, self._total_reps, self._exp_name)
        self.setup()

    @split_run.setter
    def split_run(self, value):
        '''
        Setter for split run, checks logic on selection and rebuilds submission scripts.
        :param value: bool for if run is split or not.

        :return: None
        '''
        self._split_run = value
        if self._split_run:
            # check the user has not given too many GPUS for the one replica
            if value != 1:
                raise ValueError('If you wish to run a subset of repeats per execution of TIES MD please only '
                                 'use one replica per execution and set reps_per_exec=1 in TIES.cfg')
            if len(self.devices) > 1:
                raise ValueError('1 replica per execution has been specified in TIES.cfg but multiple CUDA devices'
                                 'have been specified on the command line. Please only specify 1 device.')
            # check if we are about to run OpenMM and other instances of TIES_MD could be running we have
            # node_id set such that the output is writen to a unique location
            if self._engine == 'openmm':
                if self.node_id is None and self.run_type == 'run':
                    raise ValueError('If total_reps != reps_per_exec then the command line option --node_id'
                                     ' must be set. Please set --node_id=X where X is an integer describing which replica '
                                     'this execution of TIES_MD should be running.')
        self.sub_header, self.sub_run_line = get_header_and_run(self._engine, self.namd_version, self._split_run,
                                                                self._global_lambdas, self._total_reps, self._exp_name)

    @total_reps.setter
    def total_reps(self, value):
        '''
        Setter for total_reps, updates split_run bool.
        :param value: int for the total number of replicas.

        :return: None
        '''
        self._total_reps = value
        if self._total_reps != self._reps_per_exec:
            self.split_run = True
        else:
            self.split_run = False

    @reps_per_exec.setter
    def reps_per_exec(self, value):
        '''
        Setter for reps_per_exec, updates split_run bool.
        :param value: int for how many replicas each run of the program should execute.

        :return: None
        '''
        self._reps_per_exec = value
        if self._total_reps != self._reps_per_exec:
            self.split_run = True
        else:
            self.split_run = False

    def update_cfg(self):
        '''
        Write a TIES congig file this should be called after api changes are made.

        :return: None
        '''
        if self.engine == 'namd':
            engine = self.engine+str(self.namd_version)
        else:
            engine = self.engine
        solv_oct_box = {'cbv1': self.cell_basis_vec1[0], 'cbv2': self.cell_basis_vec1[1], 'cbv3': self.cell_basis_vec1[2],
                        'cbv4': self.cell_basis_vec2[0], 'cbv5': self.cell_basis_vec2[1], 'cbv6': self.cell_basis_vec2[2],
                        'cbv7': self.cell_basis_vec3[0], 'cbv8': self.cell_basis_vec3[1], 'cbv9': self.cell_basis_vec3[2]}
        ties_script = pkg_resources.open_text(cfg_scripts, 'TIES.cfg').read()

        ties_script = ties_script.format(engine=engine,
                                         temperature=self.temperature.in_units_of(unit.kelvin)/unit.kelvin,
                                         pressure=self.pressure.in_units_of(unit.atmospheres)/unit.atmospheres,
                                         sampling_per_window=self.sampling_per_window.in_units_of(unit.nanoseconds)/unit.nanoseconds,
                                         equili_per_window=self.equili_per_window.in_units_of(unit.nanoseconds)/unit.nanoseconds,
                                         methods=','.join(self.methods),
                                         total_reps=self.total_reps,
                                         reps_per_exec=self.reps_per_exec,
                                         elec_edges=','.join([str(x) for x in self.elec_edges]),
                                         ster_edges=','.join([str(x) for x in self.ster_edges]),
                                         global_lambdas=','.join([str(x) for x in self.global_lambdas]),
                                         cons_file='na' if self.constraint_file is None else self.constraint_file,
                                         constraint_column='na' if self.constraint_column is None else self.constraint_column,
                                         input_type=self.input_type,
                                         box_type='na' if self.box_type is None else self.box_type,
                                         edge_length='na' if self.edge_length is None else self.edge_length,
                                         **solv_oct_box)
        with open(os.path.join(self.cwd, 'TIES.cfg'), 'w') as f:
            f.write(ties_script)

    def setup(self):
        '''
        Function to setup simulations and then stop

        :return: None
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

    def build_results_dirs(self, folders):
        '''
        Helper function to build output directories.
        :param folders: List of strings for what folders to build.

        :return: None
        '''
        # If output folders do not exist make them.
        rep_dirs = range(self.total_reps)

        populated_results = list(glob.iglob(os.path.join(self.cwd, 'LAMBDA*', 'rep*', 'results', '*.npy')))
        if len(populated_results) > 0:
            raise ValueError('Results files found please check these are not needed then manually delete LAMBDA_X.XX dirs')

        exiting_lam_dirs = list(glob.iglob(os.path.join(self.cwd, 'LAMBDA*')))
        for lam_dir in exiting_lam_dirs:
            shutil.rmtree(lam_dir)

        for i in rep_dirs:
            for lam in self.lam.str_lams:
                #print('Attempting to make eq, sim and results folders for LAMBDA_{}/rep{}'.format(lam, i))
                lam_dir = 'LAMBDA_{}'.format(lam)
                for folder in folders:
                    path = os.path.join(self.cwd, lam_dir, 'rep{}'.format(i), folder)
                    Path(path).mkdir(parents=True, exist_ok=True)

    def get_options(self):
        '''
        Prints out the options the user can change and their current values.

        :return: None
        '''
        for arg in self.all_args:
            print('{}: {}'.format(arg, self.__getattribute__(arg)))

    def run(self):
        '''
        Function to run the simulations in parallel

        :return: None
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
            for lam in self.lam.str_lams:
                lam_dir = 'LAMBDA_{}'.format(lam)
                path = os.path.join(self.cwd, lam_dir, 'rep{}'.format(rep.node_id))
                if not os.path.exists(path):
                    raise ValueError('Output dir {} missing. Command line option --node_id may be set incorrectly'.format(path))

        func = partial(simulate_system, alch_sys=system, Lam=self.lam, mask=self.windows_mask,
                       cwd=self.cwd, niter=int(self.sampling_per_window/(2.0*unit.picosecond)),
                       equili_steps=int(self.equili_per_window/(2.0*unit.femtoseconds)))
        pool = Pool(processes=len(self.devices))

        tic = time.perf_counter()
        nice_print('Beginning simulations')
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

        nice_print('Finished simulations')

    def write_analysis_cfg(self):
        '''
        Function to write the input configuration files for analysis.

        :return: None
        '''

        vdw_a = ','.join(str(x) for x in self.lam.lambda_sterics_appear)
        ele_a = ','.join(str(x) for x in self.lam.lambda_electrostatics_appear)
        vdw_d = ','.join(str(x) for x in self.lam.lambda_sterics_disappear)
        ele_d = ','.join(str(x) for x in self.lam.lambda_electrostatics_disappear)

        if self.engine == 'namd':
            eng = 'NAMD'
        if self.engine == 'openmm':
            eng = 'OpenMM'

        common_cfg = '''#Temperature of simulation in units of kelvin
temperature = {0}
#Directory where any output folders are writen.
output_dir = ./analysis
#Names of thermodynamic legs in simulation, corresponds to directory names in results.
legs = {1}
#Names of engines to make analysis for (NAMD, OpenMM)
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
        #aggresivly write analysis.cfg to ensure lambdas are up to date
        with open(file_path, 'w') as f:
            f.write(common_cfg)

        file_path = os.path.join(self.cwd, '../../../exp.dat')
        if not os.path.exists(file_path):
            with open(file_path, 'w') as f:
                f.write(dummy_exp)

    def write_namd_scripts(self):
        '''
        Function to call all functions needed to write all namd input scripts

        :return: None
        '''
        self.write_namd_min()
        self.write_namd_eq()
        self.write_namd_prod()
        self.write_namd_submissions()

    def write_namd_min(self):
        '''
        Function to write eq0.conf scripts used in minimization stage

        :return: None
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
        out_name = 'sim0.conf'
        open(os.path.join(self.cwd, './replica-confs', out_name), 'w').write(min_namd_initialised)
        if not self.split_run:
            # populate and write replica script which controls replica submissions
            min_namd_uninitialised = pkg_resources.open_text(namd_sub, 'min-replicas.conf').read()
            min_namd_initialised = min_namd_uninitialised.format(reps=self.total_reps, root=self.cwd)
            out_name = 'sim0-replicas.conf'
            open(os.path.join(self.cwd, './replica-confs', out_name), 'w').write(min_namd_initialised)

    def write_namd_eq(self):
        '''
        Function to write eq1 and eq2 namd scripts used for equilibration.

        :return: None
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
            out_name = "sim{}.conf".format(i)
            open(os.path.join(self.cwd, './replica-confs', out_name), 'w').write(eq_namd_initialised)

            if not self.split_run:
                #read and write eq replica to handle replica simulations
                eq_namd_uninitialised = pkg_resources.open_text(namd_sub, 'eq-replicas.conf').read()
                eq_namd_initialised = eq_namd_uninitialised.format(reps=self.total_reps,
                                                                   prev=prev_output, current='eq{}'.format(i),
                                                                   root=self.cwd)
                open(os.path.join(self.cwd, './replica-confs', 'sim{}-replicas.conf'.format(i)), 'w').write(eq_namd_initialised)

    def write_namd_prod(self):
        '''
        Function to write sim1 namd conf script for production simulation.

        :return: None
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
        out_name = "sim3.conf"
        open(os.path.join(self.cwd, './replica-confs', out_name), 'w').write(sim_namd_initialised)

        # read and write sim replica to handle replica simulations, only if we want to use +replicas option
        if not self.split_run:
            sim_namd_uninitialised = pkg_resources.open_text(namd_sub, 'sim1-replicas.conf').read()
            sim_namd_initialised = sim_namd_uninitialised.format(reps=self.total_reps, root=self.cwd)
            open(os.path.join(self.cwd, './replica-confs', 'sim1-replicas.conf'), 'w').write(sim_namd_initialised)

    def write_namd_submissions(self):
        '''
        Function to write an example submission script of NAMD job on HPC (ARCHER2).

        :return: None
        '''
        lambs = ' '.join(self.lam.str_lams)
        if self.namd_version < 3:
            if not self.split_run:
                namd_uninitialised = pkg_resources.open_text(namd_sub, 'sub.sh').read()
                namd_initialised = namd_uninitialised.format(lambs=lambs, reps=self.total_reps, run_line=self.sub_run_line,
                                                             header=self.sub_header, root=self.cwd)
                open(os.path.join(self.cwd, 'sub.sh'), 'w').write(namd_initialised)
            else:
                namd_uninitialised = pkg_resources.open_text(namd_sub_split, 'sub.sh').read()
                namd_initialised = namd_uninitialised.format(lambs=lambs, reps=self.total_reps-1, run_line=self.sub_run_line,
                                                             header=self.sub_header, root=self.cwd)
                open(os.path.join(self.cwd, 'sub.sh'), 'w').write(namd_initialised)

        else:
            print('No NAMD3 example script currently implemented. Examples for GPU scripts can be found here '
                  'https://UCL-CCS.github.io/TIES_MD/HPC_submissions.html')

    def write_openmm_submission(self):
        '''
        Function to write an example submission script of OpenMM job on HPC (Summit).

        :return: None
        '''

        lambs = [str(x) for x in range(len(self.global_lambdas))]
        lambs = ' '.join(lambs)

        if self.split_run:
            openmm_uninitialised = pkg_resources.open_text(openmm_sub_split, 'sub.sh').read()
            openmm_initialised = openmm_uninitialised.format(lambs=lambs, reps=self.total_reps-1, run_line=self.sub_run_line,
                                                         header=self.sub_header, root=self.cwd, py_bin=sys.path[0])
            open(os.path.join(self.cwd, 'sub.sh'), 'w').write(openmm_initialised)
        else:
            # no OpenMM unified job configured for HPC
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


def nice_print(string):
    '''
    Function to print str with padding.
    :param string: A string which should be printed.

    :return: None
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
    :param exp_name: str, name of the experiment e.g. complex

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
                               ' --hint=nomultithread namd2 --tclmain sim$stage.conf $lambda $i &'

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
                               ' --hint=nomultithread namd2 +replicas {} --tclmain sim$stage-replicas.conf $lambda &'.format(reps)

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
                           ' --node_id=$i > $ties_dir/$lambda$i.out&'.format(exp_name)
        else:
            # no OpenMM unified job configured for HPC
            sub_header = None
            sub_run_line = None

    return sub_header, sub_run_line



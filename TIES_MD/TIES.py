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
#pool does not work with OpenCL
from multiprocess.pool import Pool
#ThreadPool causes inaccurate TI results if multiple threads are used
from multiprocess.pool import ThreadPool as TPool
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
from .eng_scripts import namd_sub, namd_sub_split, openmm_sub_split,\
    openmm_sub, cfg_scripts

class TIES(object):
    '''
    Class to control TIES protocol, initializes variables and calls functions to start simulation or write input scripts

    :param cwd: str, for current working directory
    :param exp_name: str, for the names of experiment i.e. complex -> complex.pdb/complex.prmtop
    :param run_type: str, flag to say if we should run dynamics or not
    :param devices: list, list of ints for which cuda devices to use
    :param rep_id: float, id denoting what replica of this simulation this execution of TIES_MD should run
    :param windows_mask: list containing ints for start and end range of windows to be run
    :param periodic: boolean determines if the simulation will be periodic
    :param lam: Lambda class, allow passing of custom lambda schedule
    :param platform: sting determines what platform OpenMM will target allowed values are ['CPU', 'CUDA', 'OpenCL', 'HIP']
    :param fast: bool, Do we want safe or fast settings for MD simulations
    :param **kwargs: dict, containing setting from config file

    '''
    def __init__(self, cwd, exp_name='complex', run_type='class', devices=None, rep_id=None, windows_mask=None,
                 periodic=True, platform='CUDA', fast=False, lam=None, **kwargs):
        nice_print('TIES')

        if run_type == 'class' and kwargs == {}:
            kwargs = read_config(os.path.join(cwd, 'TIES.cfg'))
        print('If you use this software please cite:')
        print('Wade, A.D., et al. 2022. Alchemical Free Energy Estimators and Molecular Dynamics Engines:'
              ' Accuracy, Precision, and Reproducibility. Journal of chemical theory and computation, 18(6), pp.3972-3987.\n')

        print('If you use OpenMM with this software please cite:')
        print('Chodera, J., et al. 2022. choderalab/openmmtools: 0.21.4. Zenodo. doi: 10.5281/zenodo.6625266.')
        print('Eastman, P., et al. 2017. OpenMM 7: Rapid development of high performance algorithms'
              ' for molecular dynamics. PLoS computational biology, 13(7), p.e1005659.\n')

        #check all the config file args we need are present
        args_list = ['engine', 'temperature', 'pressure', 'sampling_per_window', 'equili_per_window',
                     'methods', 'total_reps', 'split_run', 'elec_edges', 'ster_edges', 'global_lambdas', 'constraint_file',
                     'constraint_column', 'input_type', 'cell_basis_vec1', 'cell_basis_vec2', 'cell_basis_vec3']

        # check we have all required arguments
        for argument in args_list:
            if argument not in kwargs.keys():
                raise ValueError('Missing option {} in configuration file'.format(argument))

        # check we have no unexpected arguments
        for argument in kwargs.keys():
            if argument not in args_list:
                raise ValueError('Argument {} not supported for this engine or at all.'
                                 ' Please remove from the TIES.cfg.'.format(argument))

        self.all_args = args_list
        #engine must be delt with first to set namd_version which other options may need.
        api_sensitive = ['engine', 'split_run', 'elec_edges', 'ster_edges', 'global_lambdas',
                         'cell_basis_vec1', 'cell_basis_vec2', 'cell_basis_vec3']

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

        #set any attr the api needs
        self._split_run = bool(int(self._split_run))
        if self._split_run:
            self.reps_per_exec = 1
        else:
            self.reps_per_exec = self.total_reps
        self.exp_name = exp_name

        self._elec_edges = self._elec_edges.split(',')
        self._elec_edges = [float(x) for x in self._elec_edges]

        self._ster_edges = self._ster_edges.split(',')
        self._ster_edges = [float(x) for x in self._ster_edges]

        self._global_lambdas = [round(float(x), 2) for x in self._global_lambdas.split(',')]

        #set genral args
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
        self.basis_vectors = [[], [], []]

        self.total_reps = int(self.total_reps)
        self.reps_per_exec = int(self.reps_per_exec)

        if 'na' == self.constraint_file:
            self.constraint_file = None
            self.constraint_column = None

        #Deal with command line input.
        # Input for all engines is processed here if it does not pertain to the current engine its set here but not used
        self.windows = len(self._global_lambdas)
        self.windows_mask = windows_mask
        if self.windows_mask is not None:
            self.num_windows = len(range(*self.windows_mask))
        else:
            self.num_windows = self.windows

        if devices is None:
            self.devices = [0]
        else:
            if platform == 'OpenCL' or platform == 'CPU':
                raise ValueError('Currently devices option only configured for HIP and CUDA in TIES')
            self.devices = devices

        self.platform = platform
        self.cwd = cwd
        self.rep_id = rep_id
        self.periodic = periodic
        self.fast = fast

        #run through api logic
        for prop in api_sensitive:
            setattr(self, prop, self.__getattribute__('_'+prop))

        self.pre_run_line, self.run_line = None, None
        self.sub_run_line, self.sub_header = None, None

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
            TIES.setup(self)
            TIES.run(self)

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
    def cell_basis_vec1(self):
        """
        What is the 1st basis vector of the simulation cell

        :return: list of floats for x, y, z components of vector
        """
        return self._cell_basis_vec1

    @cell_basis_vec1.setter
    def cell_basis_vec1(self, value):
        '''
        Setter for cell_basis_vec1
        :param value: list for x, y, z of box, updates basis_vectors

        :return: None
        '''
        if isinstance(value, str):
            self._cell_basis_vec1 = [float(x) for x in value.split(',')]
        else:
            assert len(value) == 3
            self._cell_basis_vec1 = [x for x in value]
        self.basis_vectors[0] = Vec3(*self._cell_basis_vec1) * unit.angstrom

    @property
    def cell_basis_vec2(self):
        """
        What is the 2nd basis vector of the simulation cell

        :return: list of floats for x, y, z components of vector
        """
        return self._cell_basis_vec2

    @cell_basis_vec2.setter
    def cell_basis_vec2(self, value):
        '''
        Setter for cell_basis_vec2
        :param value: list for x, y, z of box, updates basis_vectors

        :return: None
        '''
        if isinstance(value, str):
            self._cell_basis_vec2 = [float(x) for x in value.split(',')]
        else:
            assert len(value) == 3
            self._cell_basis_vec2 = [x for x in value]
        self.basis_vectors[1] = Vec3(*self._cell_basis_vec2) * unit.angstrom

    @property
    def cell_basis_vec3(self):
        """
        What is the 3rd basis vector of the simulation cell

        :return: list of floats for x, y, z components of vector
        """
        return self._cell_basis_vec3

    @cell_basis_vec3.setter
    def cell_basis_vec3(self, value):
        '''
        Setter for cell_basis_vec3
        :param value: list for x, y, z of box, updates basis_vectors

        :return: None
        '''
        if isinstance(value, str):
            self._cell_basis_vec3 = [float(x) for x in value.split(',')]
        else:
            assert len(value) == 3
            self._cell_basis_vec3 = [x for x in value]
        self.basis_vectors[2] = Vec3(*self._cell_basis_vec3) * unit.angstrom

    @property
    def engine(self):
        '''
        What molecular dynamics engine is going to be used (openmm/namd2.14/namd3)

        :return: string for molecular dynamics engine.
        '''
        return self._engine

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
            if self.namd_version != 2.14 and self.namd_version != 3:
                raise ValueError('NAMD {} version is unsupported'
                                 ' chose from: namd2.14/namd3'.format(self.namd_version))
            if self.input_type != 'AMBER':
                raise ValueError('Only AMBER input supported in NAMD version of TIES MD')

        supported_engines = ['openmm', 'namd']
        if self._engine not in supported_engines:
            raise ValueError('Engine {} not supported please select from {}'.format(self._engine, supported_engines))

    @property
    def global_lambdas(self):
        '''
        The global values of lambda in the lambda schedule.

        :return: list of floats for lambda schedule
        '''
        return self._global_lambdas

    @global_lambdas.setter
    def global_lambdas(self, value):
        '''
        Setter for global_lambdas, rebuilds Lambdas class and submission scripts with updated info.
        :param value: list of floats for the value the global controlling parameter takes in each window.

        :return: None
        '''
        self._global_lambdas = [round(x, 2) for x in value]
        self.lam = Lambdas(self._elec_edges, self._ster_edges, self._global_lambdas, debug=False)
        self.windows = len(self._global_lambdas)

        if self.windows_mask is not None:
            try:
                check_mask = [self._global_lambdas[i] for i in range(*self.windows_mask)]
            except IndexError:
                print('Warning: windows mask range({}, {}) does not match {} windows. Change windows or'
                      ' mask.\n'.format(self.windows_mask[0], self.windows_mask[1], len(self._global_lambdas)))
        else:
            self.num_windows = self.windows

    @property
    def elec_edges(self):
        '''
        here in lambda schedule (0->1) should the electrostatic potentials begin, stop appearing.

        :return: list of floats for when the electrostatic potentials begin, stop appearing.
        '''
        return self._elec_edges

    @elec_edges.setter
    def elec_edges(self, value):
        '''
        Setter for elec_edges, rebuilds Lambdas class with updated schedule.
        :param value: list of floats, for where the electrostatic potentials begin, stop appearing.

        :return: None
        '''
        self._elec_edges = value
        self.lam = Lambdas(self._elec_edges, self._ster_edges, self._global_lambdas, debug=False)

    @property
    def ster_edges(self):
        '''
        Where in lambda schedule (0->1) should the Lennard_Jones potentials begin, stop appearing.

        :return: list of floats for when the Lennard_Jones potentials begin, stop appearing.
        '''
        return self._ster_edges

    @ster_edges.setter
    def ster_edges(self, value):
        '''
        Setter for ster_edges, rebuilds Lambdas class with updated schedule.
        :param value: list of floats, for where the Lennard_Jones potentials begin, stop appearing.

        :return: None
        '''
        self._ster_edges = value
        self.lam = Lambdas(self._elec_edges, self._ster_edges, self._global_lambdas, debug=False)

    @property
    def split_run(self):
        '''
        boolean that sets if each execution of TIESMD or NAMD runs all replicas or a subset.

        :return: bool for if run is split
        '''
        return self._split_run

    @split_run.setter
    def split_run(self, value):
        '''
        Setter for split run, checks logic on selection and rebuilds submission scripts.
        :param value: bool for if run is split or not.

        :return: None
        '''
        self._split_run = bool(value)
        if self._split_run:
            self.reps_per_exec = 1
        else:
            if self.rep_id is not None:
                raise ValueError('Split run is off but rep_id has a set value of {}. Unset rep_id'.format(self.rep_id))
            self.reps_per_exec = self.total_reps

    def update_cfg(self):
        '''
        Write a TIES config file this should be called after api changes are made.

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
                                         split_run=int(self.split_run),
                                         elec_edges=','.join([str(x) for x in self.elec_edges]),
                                         ster_edges=','.join([str(x) for x in self.ster_edges]),
                                         global_lambdas=','.join([str(x) for x in self.global_lambdas]),
                                         cons_file='na' if self.constraint_file is None else self.constraint_file,
                                         constraint_column='na' if self.constraint_column is None else self.constraint_column,
                                         input_type=self.input_type,
                                         **solv_oct_box)
        with open(os.path.join(self.cwd, 'TIES.cfg'), 'w') as f:
            f.write(ties_script)

    def setup(self):
        '''
        Function to setup simulations and then stop

        :return: None
        '''

        sub_header, sub_run_line = get_header_and_run(self.engine, self.namd_version, self.split_run,self.num_windows,
                                                      self.total_reps, self.exp_name, self.devices)
        if self.sub_header is None:
            self.sub_header = sub_header
        if self.run_line is None and self.pre_run_line is None:
            self.sub_run_line = sub_run_line
        else:
            self.sub_run_line = str(self.pre_run_line) + str(self.run_line)

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
        self.update_cfg()

    def build_results_dirs(self, folders):
        '''
        Helper function to build output directories.
        :param folders: List of strings for what folders to build.

        :return: None
        '''

        #What lambda dirs exsit
        exiting_lam_dirs = list(glob.iglob(os.path.join(self.cwd, 'LAMBDA*')))

        #do any contain results?
        pop_lam_dirs = []
        for lam_dir in exiting_lam_dirs:
            found_results = list(glob.iglob(os.path.join(lam_dir, 'rep*', 'results', '*.npy')))
            pop_lam_dirs.extend(found_results)
            found_results = list(glob.iglob(os.path.join(lam_dir, 'rep*', 'simulation', '*.alch')))
            pop_lam_dirs.extend(found_results)

        #Are the populated dirs in the schedule?
        warn_overwrite = False
        populated_lam = set()
        for lam_dir in pop_lam_dirs:
            lam_id = lam_dir.split('LAMBDA_')[1][0:4]
            populated_lam.update(lam_id)
            if lam_id not in self.lam.str_lams:
                raise ValueError('Results files found please check these are not needed'
                                 ' then manually delete LAMBDA_{} dirs'.format(lam_id))
            else:
                if self.windows_mask is not None:
                    if lam_id in [self.lam.str_lams[i] for i in range(*self.windows_mask)]:
                        warn_overwrite = True
                else:
                    if lam_id in [self.lam.str_lams[i] for i in range(self.windows)]:
                        warn_overwrite = True

        if warn_overwrite:
            print('Warning: you may be a about to re-run simulations which already have results.\n')

        # now delete anything outside of schedule assuming they are unpopulated dirs
        for lam_dir in exiting_lam_dirs:
            lam_id = lam_dir.split('LAMBDA_')[1][0:4]
            if lam_id not in self.lam.str_lams:
                print('Removing unpopulated LAMBDA_{} directory'.format(lam_id))
                shutil.rmtree('./LAMBDA_{}'.format(lam_id))

        #make all dirs needed for schedel (exist_ok=True)
        rep_dirs = range(self.total_reps)
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
        other_user_options = ['sub_header', 'pre_run_line', 'run_line']
        for arg in self.all_args+other_user_options:
            print('{}: {}'.format(arg, self.__getattribute__(arg)))

    def run(self):
        '''
        Function to run the simulations in parallel

        :return: None
        '''
        if self.split_run:
            if self.rep_id is None:
                raise ValueError('For a split run set --rep_id on the command line, or pass rep_id as an '
                                 'argument to the TIES() class.')

        system = AlchSys(self.cwd, self.exp_name, self.temperature, self.pressure, self.constraint_file,
                         self.constraint_column, self.methods, self.basis_vectors, self.input_type, self.absolute,
                         self.periodic, self.platform, self.fast)

        if self.split_run:
            system_ids = [System_ID(self.devices[0], self.rep_id)]
        else:
            system_ids = [System_ID(self.devices[i % len(self.devices)], i) for i in range(self.total_reps)]
            
        #check output dirs exist
        for rep in system_ids:
            for lam in self.lam.str_lams:
                lam_dir = 'LAMBDA_{}'.format(lam)
                path = os.path.join(self.cwd, lam_dir, 'rep{}'.format(rep.rep_id))
                if not os.path.exists(path):
                    raise ValueError('Output dir {} missing. --rep_id may be incorrect.'.format(path))

        #make mask so all windows are run if user does not pass mask.
        if self.windows_mask is None:
            self.windows_mask = [0, self.windows]

        nice_print('Beginning simulations')
        func = partial(simulate_system, alch_sys=system, Lam=self.lam, mask=self.windows_mask,
                       cwd=self.cwd, niter=int(self.sampling_per_window/(2.0*unit.picosecond)),
                       equili_steps=int(self.equili_per_window/(2.0*unit.femtoseconds)))
        if self.platform != 'OpenCL':
            pool = Pool(processes=len(self.devices))
        else:
            pool = TPool(processes=len(self.devices))

        tic = time.perf_counter()
        simulated = pool.map(func, system_ids) ### MAIN MD
        toc = time.perf_counter()
        total_simulation_time = (toc - tic)

        # Kill extra threads
        pool.close()
        pool.join()
        pool.terminate()

        total_sampling = (((self.sampling_per_window + self.equili_per_window)
                           * self.num_windows * self.reps_per_exec) / len(set(self.devices))) / unit.nanosecond

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

        file_path = os.path.join(self.cwd, '../../../analysis.cfg')
        #aggresivly write analysis.cfg to ensure lambdas are up to date
        with open(file_path, 'w') as f:
            f.write(common_cfg)

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
        out_name = 'run0.conf'
        open(os.path.join(self.cwd, './replica-confs', out_name), 'w').write(min_namd_initialised)
        if not self.split_run:
            # populate and write replica script which controls replica submissions
            min_namd_uninitialised = pkg_resources.open_text(namd_sub, 'min-replicas.conf').read()
            min_namd_initialised = min_namd_uninitialised.format(reps=self.total_reps, root=self.cwd)
            out_name = 'run0-replicas.conf'
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

                # Older TIES protocol uses Berendsen (no longer supported or tested)
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
            out_name = "run{}.conf".format(i)
            open(os.path.join(self.cwd, './replica-confs', out_name), 'w').write(eq_namd_initialised)

            if not self.split_run:
                #read and write eq replica to handle replica simulations
                eq_namd_uninitialised = pkg_resources.open_text(namd_sub, 'eq-replicas.conf').read()
                eq_namd_initialised = eq_namd_uninitialised.format(reps=self.total_reps,
                                                                   prev=prev_output, current=i,
                                                                   root=self.cwd)
                open(os.path.join(self.cwd, './replica-confs', 'run{}-replicas.conf'.format(i)), 'w').write(eq_namd_initialised)

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
        out_name = "run3.conf"
        open(os.path.join(self.cwd, './replica-confs', out_name), 'w').write(sim_namd_initialised)

        # read and write sim replica to handle replica simulations, only if we want to use +replicas option
        if not self.split_run:
            sim_namd_uninitialised = pkg_resources.open_text(namd_sub, 'sim1-replicas.conf').read()
            sim_namd_initialised = sim_namd_uninitialised.format(reps=self.total_reps, root=self.cwd)
            open(os.path.join(self.cwd, './replica-confs', 'run3-replicas.conf'), 'w').write(sim_namd_initialised)

    def write_namd_submissions(self):
        '''
        Function to write an example submission script of NAMD job on HPC (ARCHER2).

        :return: None
        '''
        if self.windows_mask is not None:
            lambs = [self.lam.str_lams[i] for i in range(*self.windows_mask)]
        else:
            lambs = self.lam.str_lams
        lambs = ' '.join(lambs)

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

        if self.windows_mask is not None:
            lambs = [str(i) for i in range(*self.windows_mask)]
        else:
            lambs = [str(i) for i in range(self.windows)]
        lambs = ' '.join(lambs)

        if self.split_run:
            openmm_uninitialised = pkg_resources.open_text(openmm_sub_split, 'sub.sh').read()
            openmm_initialised = openmm_uninitialised.format(lambs=lambs, reps=self.total_reps-1, run_line=self.sub_run_line,
                                                         header=self.sub_header, root=self.cwd, py_bin=sys.path[0])
            open(os.path.join(self.cwd, 'sub.sh'), 'w').write(openmm_initialised)
        else:
            openmm_uninitialised = pkg_resources.open_text(openmm_sub, 'sub.sh').read()
            openmm_initialised = openmm_uninitialised.format(lambs=lambs, run_line=self.sub_run_line,
                                                             header=self.sub_header, root=self.cwd, py_bin=sys.path[0])
            open(os.path.join(self.cwd, 'sub.sh'), 'w').write(openmm_initialised)

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


def get_header_and_run(engine, namd_version, split_run, num_windows, reps, exp_name, devices):
    '''
    Function to prep submission file. Number of windows and replicas are inspected to make best guess at
    number of nodes and CPUS/GPUS

    :param engine: str, What engine are we using [namd, openmm]
    :param namd_version: float, What version of namd are we using if any
    :param split_run: bool, Should each run line run all or one replica
    :param num_windows: int for number of windows running
    :param reps: int, number of replica simulations
    :param exp_name: str, name of the experiment e.g. complex
    :param devices: list of ints for which gpus to target

    :return: str, str for header and run line of HPC sub script
    '''

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
                               ' --hint=nomultithread namd2 --tclmain run$stage.conf $lambda $i &'

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

                sub_run_line = 'srun -N $nodes_per_namd -n $cpus_per_namd --distribution=block:block --hint=nomultithread' \
                               ' namd2 +replicas {} --tclmain run$stage-replicas.conf $lambda &'.format(reps)

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
            sub_run_line = 'jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg' \
                           ' --exp_name={} --windows_mask=$lambda,$(expr $lambda + 1)' \
                           ' --rep_id=$i > $ties_dir/$lambda_$i.out&'.format(exp_name)
        else:
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
#BSUB -e eLIGPAIR.%J""".format(int(np.ceil(num_jobs / gpus_per_node)))
            sub_run_line = 'jsrun --smpiargs="off" -n 1 -a 1 -c {} -g {} -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg' \
                           ' --exp_name={} --windows_mask=$lambda,$(expr $lambda + 1) --devices={}' \
                           ' > $ties_dir/$lambda_$i.out&'.format(len(devices), len(devices), exp_name, ','.join([str(x) for x in devices]))

    return sub_header, sub_run_line

def read_config(config_file):
    '''
    Function to read config file from disk
    :param config_file: str pointing to TIES.cfg file
    :return: dict, containing all config file args
    '''

    args_dict = {}
    with open(config_file) as file:
        for line in file:
            if line[0] != '#' and line[0] != ' ' and line[0] != '\n':
                data = line.rstrip('\n').split('=')
                if len(data) > 2:
                    raise ValueError('Failed to parse line: {}'.format(line))
                # Remove spaces
                data = [s.replace(" ", "") for s in data]
                #remove tabs
                data = [s.replace("\t", "") for s in data]
                args_dict[data[0]] = data[1]
    return args_dict



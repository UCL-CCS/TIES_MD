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

import os
import glob
import collections
import numpy as np

from ..methods.TI import TI_Analysis
from ..methods.FEP import MBAR_Analysis

class Lambdas():
    '''
    Class: holds the information about a lambda schedule in an easily queryable format

    :param vdw_a: list, Contains floats for values of VDW appearing lambdas
    :param vdw_d: list, Contains floats for values of VDW disappearing lambdas
    :param ele_a: list, Contains floats for values of ELEC appearing lambdas
    :param ele_d: list, Contains floats for values of ELEC disappearing lambdas
    '''

    def __init__(self, vdw_a, vdw_d, ele_a, ele_d):

        self.lambda_sterics_appear = [float(x) for x in vdw_a]
        self.lambda_electrostatics_appear = [float(x) for x in ele_a]
        self.lambda_sterics_disappear = [float(x) for x in vdw_d]
        self.lambda_electrostatics_disappear = [float(x) for x in ele_d]

        assert (len(self.lambda_sterics_appear) == len(self.lambda_electrostatics_appear))
        assert (len(self.lambda_sterics_disappear) == len(self.lambda_electrostatics_appear))
        assert (len(self.lambda_sterics_disappear) == len(self.lambda_electrostatics_disappear))

        self.schedule = []
        for i, (su, sd, eu, ed) in enumerate(zip(self.lambda_sterics_appear, self.lambda_sterics_disappear,
                                                 self.lambda_electrostatics_appear,
                                                 self.lambda_electrostatics_disappear)):
            param_vals = {'lambda_sterics_appear': su, 'lambda_sterics_disappear': sd,
                          'lambda_electrostatics_appear': eu, 'lambda_electrostatics_disappear': ed}
            self.schedule.append(param_vals)

    def update_attrs_from_schedule(self):
        '''
        helper function to update the values of self.lambda_sterics_appear etc if the self.schedule is changed
        '''
        self.lambda_sterics_appear = [x['lambda_sterics_appear'] for x in self.schedule]
        self.lambda_electrostatics_appear = [x['lambda_electrostatics_appear'] for x in self.schedule]
        self.lambda_sterics_disappear = [x['lambda_sterics_disappear'] for x in self.schedule]
        self.lambda_electrostatics_disappear = [x['lambda_electrostatics_disappear'] for x in self.schedule]

class OpenMM(object):
    '''
    Class to perform TIES analysis on OpenMM_TIES results

    :param method: str, 'TI' or 'FEP'
    :param output: str, pointing to base dir of where output will be writen
    :param win_mask: list of ints, what windows if any to remove from analysis
    :param distributions: bool, Do we want to calculate the dG for each rep individually
    :param rep_convg: list of ints, what intermediate number of reps do you wish to inspect convergence for
    :param sampling_convg: list of ints, what intermediate amount of sampling do
     you wish to inspect convergence for
    :param vdw_a: list of floats, describes lambda schedule for vdw appear
    :param vdw_d: list of floats, describes lambda schedule for vdw disappear
    :param ele_a: list of floats, describes lambda schedule for elec appear
    :param ele_d: list of floats, describes lambda schedule for elec disappear
    :param fep_combine_reps: bool: 1 or 0 an option to combine fep replicas into one time series
    '''
    def __init__(self, method, output, win_mask, distributions, rep_convg, sampling_convg,
                 vdw_a, vdw_d, ele_a, ele_d, fep_combine_reps):
        self.name = 'OpenMM'
        self.method = method
        self.openmm_lambs = Lambdas(vdw_a, vdw_d, ele_a, ele_d)
        self.output = output
        self.fep_combine_reps = bool(int(fep_combine_reps))
        self.win_mask = win_mask
        self.distributions = distributions
        self.rep_convg = rep_convg
        self.sampling_convg = sampling_convg

    def run_analysis(self, data_root, temp, prot, lig, leg):
        '''
        Function to run the analysis for each method allowed for this engine.

        :param data_root: str, file path point to base dir for results files
        :param temp: float for temperature in units of kelvin
        :param prot: str, name of dir for protein
        :param lig:  str, name of dir for ligand
        :param leg: str, name of dir for thermo leg
        :return: list of floats, [dg, stdev(dg)]
        '''

        analysis_dir = os.path.join(self.output, self.name, self.method, prot, lig, leg)
        if self.method == 'FEP':
            try:
                data = self.collate_data(data_root, prot, lig, leg)
            except ValueError:
                raise ValueError('If you are attempting to add extra windows to FEP(MBAR) this does not work try TI only.')
            method_run = MBAR_Analysis(data, temp, self.openmm_lambs, analysis_dir)
        elif self.method == 'TI':
            data = self.collate_data(data_root, prot, lig, leg)
            method_run = TI_Analysis(data, self.openmm_lambs, analysis_dir)
        else:
            raise ValueError('Unknown method {}'.format(self.method))

        #provide option to combine all decorrelated timeseries into one long psuedo trajectory
        if self.method == 'FEP' and self.fep_combine_reps == True:
            print('Running MBAR analysis with all replica time series combined.')
            print('Errors will be estimated by PYMBAR')
            if self.distributions:
                print('If combining FEP data distributions cant be calculated ignoring distributions=True in analysis.cfg')
            result = method_run.analysis(mask_windows=self.win_mask)
            #convert stdev given by pymbar to SEM
            result = [result[0], result[1]/np.sqrt(data.shape[0])]
        elif self.method == 'FEP':
            result = method_run.replica_analysis(self.distributions, self.rep_convg, self.sampling_convg, self.win_mask)

        else:
            #TI analysis
            result = method_run.analysis(self.distributions, self.rep_convg, self.sampling_convg, self.win_mask)

        return result

    def collate_data(self, data_root, prot, lig, leg):
        '''
        Function to iterate over replica and window dirs reading OpenMM outputs and building numpy array of potentials

        :param data_root: str, file path point to base dir for results files
        :param prot: str, name of dir for protein
        :param lig:  str, name of dir for ligand
        :param leg: str, name of dir for thermo leg

        :return: np.array() containing all the data in the result files concatenated
        '''

        results_dir_path = os.path.join(data_root, prot, lig, leg)

        if self.method == 'TI':
            TI_results_files = os.path.join(results_dir_path,  'LAMBDA_*', 'rep*', 'results',  '*TI.npy')
            result_files = list(glob.iglob(TI_results_files))
        elif self.method == 'FEP':
            FEP_results_files = os.path.join(results_dir_path, 'LAMBDA_*', 'rep*', 'results', '*FEP.npy')
            result_files = list(glob.iglob(FEP_results_files))
        else:
            raise ValueError('Unknown method specified, {}, choose from TI/FEP'.format(self.method))

        if len(result_files) == 0:
            raise ValueError('{} in methods but no results files found'.format(self.method))

        # Sort by order of windows
        result_files.sort(key=get_replica)
        result_files.sort(key=get_window)

        # Use ordered dict to preserve windows order
        all_data = collections.OrderedDict()
        for file in result_files:
            window = get_window(file)
            loaded = np.load(file, allow_pickle=True)
            loaded = np.array([loaded])
            if window not in all_data:
                all_data[window] = loaded
            else:
                # appending all repeats to a new axis
                all_data[window] = np.vstack([all_data[window], loaded])

        # appending all windows together to make final array
        concat_windows = np.stack([x for x in all_data.values()], axis=0)
        # resuffle axis to be in order reps, windows, lambda_dimensions, iterations
        concat_windows = np.transpose(concat_windows, (1, 0, 2, 3))

        return concat_windows

def get_window(string):
    '''
    Helper function to sort directory paths by specific index in file name

    :param string: File path to results file

    :return: float for the window value i.e. LAMBDA_0.00 return 0.00
    '''
    path = os.path.normpath(string)
    return float(path.split(os.sep)[-4].split('_')[1])

def get_replica(string):
    '''
    Helper function to sort directory paths by specific index in file name

    :param string: File path to results file

    :return: int for the replica id i.e. rep0 return 0
    '''
    path = os.path.normpath(string)
    return int(path.split(os.sep)[-3].split('rep')[1])

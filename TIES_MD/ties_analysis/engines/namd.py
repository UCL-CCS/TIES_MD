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

from .openmm import Lambdas
from ..methods.TI import TI_Analysis

class NAMD(object):
    '''
    Class to perform TIES analysis on NAMD results

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
    :param namd_version: float, for which version of NAMD generated alch files
    '''
    def __init__(self, method, output, win_mask, distributions, rep_convg, sampling_convg,
                 vdw_a, vdw_d, ele_a, ele_d, namd_version):
        self.namd_ver = float(namd_version)

        if self.namd_ver < 3:
            self.name = 'NAMD2'
        else:
            self.name = 'NAMD3'
        self.method = method
        self.namd_lambs = Lambdas(vdw_a, vdw_d, ele_a, ele_d)
        self.output = output
        self.win_mask = win_mask
        self.distributions = distributions
        self.rep_convg = rep_convg
        self.sampling_convg = sampling_convg

    def run_analysis(self,  data_root, temp, prot, lig, leg):
        '''
        Function to run the analysis for each method allowed for this engine.

        :param data_root: str, file path point to base dir for results files
        :param temp: float for temperature in units of kelvin (not used in NAMD as FEP not implemented)
        :param prot: str, name of dir for protein
        :param lig:  str, name of dir for ligand
        :param leg: str, name of dir for thermo leg

        :return: list of floats, [dg, stdev(dg)]
        '''

        if self.method == 'FEP':
            raise NotImplementedError('FEP not supported in NAMD analysis currently.')

        data = self.collate_data(data_root, prot, lig, leg)
        analysis_dir = os.path.join(self.output, self.name, self.method, prot, lig, leg)

        method_run = TI_Analysis(data, self.namd_lambs, analysis_dir)
        result = method_run.analysis(self.distributions, self.rep_convg, self.sampling_convg, self.win_mask)

        return result

    def collate_data(self, data_root, prot, lig, leg):
        '''
        Function to iterate over replica and window dirs reading NAMD alch file and building numpy array of potentials

        :param data_root: str, file path point to base dir for results files
        :param prot: str, name of dir for protein
        :param lig:  str, name of dir for ligand
        :param leg: str, name of dir for thermo leg

        :return: numpy.array for all potential collected
        '''

        results_dir_path = os.path.join(data_root, prot, lig, leg)

        result_files = os.path.join(results_dir_path, 'LAMBDA_*', 'rep*', 'simulation', 'sim1.alch')
        result_files = list(glob.iglob(result_files))

        if len(result_files) == 0:
            raise ValueError('{} in methods but no results files found'.format(self.method))

        # Sort by order of replicas then windows
        result_files.sort(key=get_replica)
        result_files.sort(key=get_window)

        iterations = get_iter(result_files[0])

        #print('Processing files...')
        #for file in result_files:
        #    print(file)

        # Use ordered dict to preserve windows order
        all_data = collections.OrderedDict()
        for file in result_files:
            window = get_window(file)
            data = read_alch_file(file, self.namd_ver, iterations)
            data = np.array([data])
            if window not in all_data:
                all_data[window] = data
            else:
                #stack reps of same window together
                all_data[window] = np.vstack([all_data[window], data])

        # appending all windows together to make final array
        concat_windows = np.stack([x for x in all_data.values()], axis=0)
        #resuffle axis to be in order reps, windows, lambda_dimensions, iterations
        concat_windows = np.transpose(concat_windows, (1, 0, 2, 3))

        return concat_windows

def read_alch_file(file_path, namd_ver, iterations):
    '''
    Function for reading different NAMD ver. alch files

    :param file_path: str, location of namd alch file
    :param namd_ver: float, new or old used to specify what format of namd alch file we are looking at (old <= 2.12)
    :param iterations: int, Number sample in alch file

    :return: numpy array, contains potentials from one namd alch file
    '''
    # final data has order sterics appear/dis elec appear/dis to match openmm
    data = np.zeros([4, iterations])
    with open(file_path) as f:
        count = 0
        for line in f:
            if line[0:2] == 'TI':
                split_line = line.split()
                if namd_ver > 2.12:
                    data[:, count] = [float(split_line[6]), float(split_line[12]),
                                      float(split_line[4]), float(split_line[10])]
                elif namd_ver <= 2.12:
                    data[:, count] = [float(split_line[4]), float(split_line[8]),
                                      float(split_line[2]), float(split_line[6])]
                else:
                    raise ValueError('Unknown NAMD ver. {}'.format(namd_ver))

                count += 1
    if count != iterations:
        print('WARNING: {} terminated early only found {}/{} iterations.'.format(file_path, count, iterations))
    return data

def get_iter(file_loc):
    '''
    Function to get the number of gradient samples in an NAMD alch file

    :param file_loc: file path to alch file

    :return: int for the number of gradient samples
    '''
    iterations = 0
    with open(file_loc) as f:
        for line in f:
            if line[0:2] == 'TI':
                iterations += 1
    return iterations


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

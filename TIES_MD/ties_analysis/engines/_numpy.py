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
from ..methods.FEP import MBAR_Analysis


class Numpy(object):
    '''
    Class to perform TIES analysis on numpy array

    :param method: str, 'TI' or 'FEP'
    :param output: str, pointing to base dir of where output will be writen
    :param win_mask: list of ints, what windows if any to remove from analysis
    :param distributions bool, Do we want to calculate the dG for each rep individually
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

        self.name = 'Numpy'
        self.method = method
        self.lambs = Lambdas(vdw_a, vdw_d, ele_a, ele_d)
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

        data = np.load(os.path.join(data_root, prot, lig, leg, 'results.npy'))
        analysis_dir = os.path.join(self.output, self.name, self.method, prot, lig, leg)

        if self.method == 'FEP':
            method_run = MBAR_Analysis(data, temp, self.lambs, analysis_dir)
        elif self.method == 'TI':
            method_run = TI_Analysis(data, self.lambs, analysis_dir)
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

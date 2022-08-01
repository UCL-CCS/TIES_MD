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

import numpy as np
import copy
import os
from pathlib import Path

from pymbar import MBAR, timeseries

from ..methods.TI import compute_bs_error

skip_graphs = False
try:
    import seaborn as sns
    from matplotlib import pyplot as plt
    from matplotlib import colors

    sns.set_style('whitegrid')
    sns.set_context('paper', font_scale=2.5)

except ModuleNotFoundError:
    print('matplotlib/seaborn not found skipping graphs')
    skip_graphs = True

#GLOBAL CONSTANTS
kb = 0.0019872066135803576 #unit kcal/(mole*kelvin)

class MBAR_Analysis():
    '''
    Class for MBAR analysis

    :param MBARs: numpy array for all results
    :param temp: float for temperature in units of kelvin
    :param lambdas: Lambda class containing schedule
    :param distribution: Boolean, if True then dGs will not be averaged and a distribution of results is returned
    :param analysis_dir: string file path for where to save analysis output
    :param decorrelate: bool, do we want to decorrelate data before processing.
    '''

    def __init__(self, MBARs, temp, lambdas, analysis_dir, decorrelate=False):

        self.data = MBARs

        self.KBT = kb*temp #unit kcal/mol
        self.shape = list(self.data.shape)
        self.nstates = self.shape[1]
        print('Data shape is {} repeats, {} states_i, {} states_j, {} iterations\n'.format(*self.shape))

        self.lambdas = lambdas

        self.analysis_dir = analysis_dir
        if analysis_dir is not None:
            print('Attempting to make analysis folder {0}'.format(self.analysis_dir))
            Path(self.analysis_dir).mkdir(parents=True, exist_ok=True)

        self.data, self.N_k, self.rep_data, self.rep_N_k = MBAR_Analysis.decorrelate_data(self, decorrelate)

    def decorrelate_data(self, decorrelate):
        '''
        Decorrolate time series data.

        :param decorrelate: boolean, True if we want decorrelated data else False

        :return: turple, (np matrix containing decorrelated data, list of ints for end of decorrelated data in matrix)
        '''

        #load data if avaliable
        if self.analysis_dir is not None:
            decorr_data_sav = os.path.join(self.analysis_dir, 'fep_decorr_data.npy')
            index_sav = os.path.join(self.analysis_dir, 'fep_N_k.npy')

            rep_decorr_data_sav = os.path.join(self.analysis_dir, 'fep_rep_decorr_data.npy')
            rep_index_sav = os.path.join(self.analysis_dir, 'fep_rep_N_k.npy')
            
            """
            if os.path.exists(decorr_data_sav) and os.path.exists(index_sav)\
                    and os.path.exists(rep_decorr_data_sav) and os.path.exists(rep_index_sav):

                print('Loading decorrelated data from disk')
                decorr_data = np.load(decorr_data_sav)
                N_k = np.load(index_sav)

                replicas_deccor = np.load(rep_decorr_data_sav)
                replicas_Nk = np.load(rep_index_sav)
                return decorr_data, N_k, replicas_deccor, replicas_Nk
            """

        print('Decorrelating data...')
        iter_per_rep = self.shape[3]
        num_repeats = self.shape[0]
        tot_iter = iter_per_rep * num_repeats

        # decorr data matrix is square for numpy convenience and padded with zeros.
        # N_k allows us to know were the data ends and the padding begins in the matrix.
        N_k = np.zeros([self.nstates], np.int32)
        decorr_data = np.zeros([self.nstates, self.nstates, tot_iter])

        # to compute ensemble stdev we also would like to save the result for each replica separately
        replicas_deccor = []
        replicas_Nk = []
        blank_decorr_data = np.zeros([self.nstates, self.nstates, iter_per_rep])
        blank_Nk = np.zeros([self.nstates], np.int32)

        for i, u_kln in enumerate(self.data):
            if decorrelate:
                rep_deccor_data = copy.deepcopy(blank_decorr_data)
                rep_Nk = copy.deepcopy(blank_Nk)
                for k in range(self.nstates):
                        [nequil, g, Neff_max] = timeseries.detect_equilibration(u_kln[k, k, :])
                        sub_idx = timeseries.subsample_correlated_data(u_kln[k, k, :], g=g)
                        decorr_data[k, :, 0 + N_k[k]:N_k[k] + len(sub_idx)] = u_kln[k, :, sub_idx].T
                        rep_deccor_data[k, :, 0:len(sub_idx)] = u_kln[k, :, sub_idx].T
                        N_k[k] += len(sub_idx)
                        rep_Nk[k] = len(sub_idx)
            else:
                rep_deccor_data = u_kln
                rep_Nk = [self.shape[-1] for i in range(self.nstates)]
                N_k = [x+self.shape[-1] for x in N_k]

            replicas_deccor.append(rep_deccor_data)
            replicas_Nk.append(rep_Nk)

        if decorrelate:
            print('Number of decorrelated samples extracted per state: {}'.format(N_k))
        else:
            print('Decorrelatation was skipped using max iters per state: {}'.format(N_k))

        #if self.analysis_dir is not None:
        #    np.save(decorr_data_sav, decorr_data)
        #    np.save(index_sav, N_k)

        #    np.save(rep_decorr_data_sav, replicas_deccor)
        #    np.save(rep_index_sav, replicas_Nk)

        return decorr_data, N_k, replicas_deccor, replicas_Nk

    def replica_analysis(self, distributions=False, rep_convg=None, sampling_convg=None, mask_windows=None):
        '''
        Function to make analysis of result from MBAR considering each trajectory as one replica

        :param distributions: bool, Do we want to calculate the dG for each rep individually
        :param rep_convg: list of ints, what number of reps do we want results for.
        :param sampling_convg: list of ints, what number of samples do we want results for.
        :param mask_windows: list of ints, can be used to specify what windows to remove.

        :return: containing average of bootstrapped dG and SEM
        '''

        if sampling_convg is not None:
            sampling_free_energy = []
            for sampling in sampling_convg:
                if sampling > self.shape[-1]:
                    raise ValueError('Requested convergence info for too large iteration count: {}.'
                                     ' Max number of iterations is {}'.format(sampling, self.shape[-1]))
                replica_results = []
                for d, n in zip(self.rep_data, self.rep_N_k):
                    tmp_nk = [x if x < sampling else sampling for x in n]
                    mbar_res = self.analysis(u_kln=d[:, :, 0:sampling], N_k=tmp_nk, mask_windows=mask_windows)
                    replica_results.append(mbar_res)
                dg = [x[0] for x in replica_results]
                bs_res = compute_bs_error(dg)
                sampling_free_energy.append([bs_res[0], np.sqrt(bs_res[1])])
            print('Convergence with number of samples:')
            print(sampling_convg)
            print(sampling_free_energy)
            print('')

        if rep_convg is not None:
            rep_free_energy = []
            for rep in rep_convg:
                replica_results = []
                for d, n in zip(self.rep_data[0:rep], self.rep_N_k[0:rep]):
                    mbar_res = self.analysis(u_kln=d, N_k=n, mask_windows=mask_windows)
                    replica_results.append(mbar_res)
                dg = [x[0] for x in replica_results]
                bs_res = compute_bs_error(dg)
                rep_free_energy.append([bs_res[0], np.sqrt(bs_res[1])])
            print('Convergence with number of reps:')
            print(rep_convg)
            print(rep_free_energy)
            print('')

        replica_results = []
        for i, (d, n) in enumerate(zip(self.rep_data, self.rep_N_k)):
            mbar_res = self.analysis(u_kln=d, N_k=n, mask_windows=mask_windows, rep_id=i)
            replica_results.append(mbar_res)

        dg = [x[0] for x in replica_results]

        if distributions:
            print('dG for each replica:')
            print(dg)
            mbar_err = np.average([x[1] / np.sqrt(self.shape[0]) for x in replica_results])
            print('With average MBAR error: {}'.format(mbar_err))
            print('')

        bs_res = compute_bs_error(dg)
        return [bs_res[0], np.sqrt(bs_res[1])]

    def analysis(self, u_kln=None, N_k=None, mask_windows=None, rep_id=None):
        '''
        Process a matrix of potentials passes to this function as u_kln or process self.data which is matrix of all replicas

        :param u_kln: numpy array matrix of potentials for one replica
        :param N_k: list of ints, specifies which entry in u_kln should be read up to for each window
        :param mask_windows: list of ints, can be used to specify what windows to remove from u_kln

        :return: list, containing dG and associated error calculated by MBAR
        '''

        if u_kln is None:
            u_kln = self.data
        if N_k is None:
            N_k = self.N_k

        # remove windows
        if mask_windows is not None:
            mask_windows.sort()
            mask_windows.reverse()
            print('Removing windows {}'.format(mask_windows))
            for win in mask_windows:
                u_kln = np.delete(u_kln, win, axis=0)
                u_kln = np.delete(u_kln, win, axis=1)
                N_k = np.delete(N_k, win, axis=0)

        # Compute free energy differences and statistical uncertainties
        mbar = MBAR(u_kln, N_k)
        mbar_results = mbar.compute_free_energy_differences(return_theta=True)

        # print("Number of uncorrelated samples per state: {}".format(N_k))
        result = ([mbar_results['Delta_f'][0, len(N_k) - 1]*self.KBT,
                   mbar_results['dDelta_f'][0, len(N_k) - 1]*self.KBT]) #units kilocalorie per mol

        # Save data to analysis dir
        if self.analysis_dir is not None:
            np.save(os.path.join(self.analysis_dir, 'dG_by_state.npy'), mbar_results['Delta_f'])
            #np.save(os.path.join(self.analysis_dir, 'sigma_dG_by_state.npy'), mbar_results['dDelta_f'])
            if rep_id is not None:
                over_lap_res = mbar.compute_overlap()
                np.save(os.path.join(self.analysis_dir, 'overlap{}.npy'.format(rep_id)), over_lap_res['matrix'])
                if not skip_graphs:
                    self.plot_overlap_mat(over_lap_res['matrix'], rep_id)

        return result

    def plot_overlap_mat(self, mat, rep_id):
        '''
        Make a plot of the overlap matrix for this simulation

        :param mat: numpy array, overlap matrix to make a plot for
        :param rep_id: int, what replica are we looking at
        '''
        print('Plotting overlap matrix for replica {}...'.format(rep_id))
        plt.rcParams['figure.figsize'] = self.shape[1], self.shape[1]
        plt.rcParams["font.size"] = 22
        plt.rcParams["font.serif"] = "Computer Modern Roman"
        plt.rcParams["mathtext.fontset"] = "cm"

        cmap = colors.ListedColormap(['#FBE8EB', '#88CCEE', '#78C592', '#117733'])
        bounds = [0.0, 0.025, 0.1, 0.3, 0.8]
        norm = colors.BoundaryNorm(bounds, cmap.N, clip=False)
        cbar_kws = dict(ticks=[.025, .1, .3, 0.8])
        ax = sns.heatmap(mat, annot=True, fmt='.2f', linewidths=.3, annot_kws={"size": 14},
                         square=True, robust=True, cmap=cmap, norm=norm, vmin=0, vmax=1, cbar_kws=cbar_kws)
        ax.set_xlabel(r'$\lambda$ index')
        ax.set_ylabel(r'$\lambda$ index')
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.set_xlabel(r'$\lambda$ index')
        ax.set_ylabel(r'$\lambda$ index')

        # Plot
        plt.tight_layout()
        plt.savefig(os.path.join(self.analysis_dir, 'overlap_for_rep{}.png'.format(rep_id)))
        plt.close('all')

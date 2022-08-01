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

import copy
import numpy as np
import os
from pathlib import Path
from sklearn.utils import resample

skip_graphs = False
try:
    from matplotlib import pyplot as plt
    from matplotlib import colors

except ModuleNotFoundError:
    print('matplotlib not found skipping graphs')
    skip_graphs = True

class TI_Analysis(object):
    '''
    Class for thermodynamic integration analysis

    :param grads: numpy array for all results
    :param lambdas: Lambda class containing schedule
    :param distribution: Boolean, if True then dGs will not be averaged and a distribution of results is returned
    :param analysis_dir: string, pointing to where we want analysis saved
    '''
    def __init__(self, grads, lambdas, analysis_dir):

        self.data = grads

        self.lambdas = lambdas
        self.shape = list(self.data.shape)
        self.nstates = self.shape[1]

        #directory where any TI specific output could be saved, not curretly used
        self.analysis_dir = analysis_dir
        if analysis_dir is not None:
            print('Attempting to make analysis folder {0}'.format(self.analysis_dir))
            Path(self.analysis_dir).mkdir(parents=True, exist_ok=True)

        print('Data shape is {} repeats, {} states, {} lambda dimensions, {} iterations\n'.format(*self.shape))
        if not skip_graphs:
            self.plot_du_by_dl()

    def analysis(self, distributions=False, rep_convg=None, sampling_convg=None, mask_windows=None):
        '''
        Perform TI analysis

        :param distributions: bool, Do we want to calculate the dG for each rep individually
        :param rep_convg: list of ints, what number of reps do we want results for.
        :param sampling_convg: list of ints, what number of samples do we want results for.
        :param mask_windows: list of ints, can be used to specify what windows to remove.

        :return: list, dG calculated by TI and associated standard deviation
        '''

        data = copy.deepcopy(self.data)
        lambdas = copy.deepcopy(self.lambdas)

        # remove windows if requested
        if mask_windows is not None:
            mask_windows.sort()
            mask_windows.reverse()
            print('Removing windows {}'.format(mask_windows))
            for win in mask_windows:
                data = np.delete(data, win, axis=1)
                del lambdas.schedule[win]
        lambdas.update_attrs_from_schedule()

        #iterate over sampling and compute dg if requested
        if sampling_convg is not None:
            sampling_free_energy = []
            for sampling in sampling_convg:
                if sampling > self.shape[-1]:
                    raise ValueError('Requested convergence info for too large iteration count: {}.'
                                     ' Max number of iterations is {}'.format(sampling, self.shape[-1]))
                avg_over_iterations = np.average(data[:, :, :, 0:sampling], axis=3)
                free_energy, variance = self.intergrate(avg_over_iterations.transpose(), lambdas)
                result = [sum(free_energy.values()), np.sqrt(sum(variance.values()))]
                sampling_free_energy.append(result)
            print('Convergence with number of samples:')
            print(sampling_convg)
            print(sampling_free_energy)
            print('')

        avg_over_iterations = np.average(data, axis=3)
        #compute dgs seperatly for each replica if requested
        if distributions:
            dist_free_energy = []
            for rep in avg_over_iterations:
                free_energy, _ = self.intergrate(np.array([rep]).transpose(), lambdas)
                dist_free_energy.append(sum(free_energy.values()))
            print('dG for each replica:')
            print(dist_free_energy)
            print('')

        #iterate overs chuncks of replicas and compute dg if requested
        if rep_convg is not None:
            rep_free_energy = []
            for rep in rep_convg:
                free_energy, variance = self.intergrate(avg_over_iterations[0:rep].transpose(), lambdas)
                result = [sum(free_energy.values()), np.sqrt(sum(variance.values()))]
                rep_free_energy.append(result)
            print('Convergence with number of reps:')
            print(rep_convg)
            print(rep_free_energy)
            print('')

        #compute dg for all reps and sampling this is the result that is returned
        free_energy, variance = self.intergrate(avg_over_iterations.transpose(), lambdas)

        #print('Free energy breakdown:')
        #print(free_energy)
        #print('Variance breakdown:')
        #print(variance)
        #print('')

        return [sum(free_energy.values()), np.sqrt(sum(variance.values()))]

    def intergrate(self, data, lambdas):
        '''
        Function to perform numerical integration.

        :param data: numpy array of averaged gradients
        :param lambdas: class, contains information about lambda schedule

        :return: turple of dicts, {lambda_parameters: free energy}, {lambda_parameters: free energy variance} for
         each lambda parameter
        '''

        lambda_names = lambdas.schedule[0].keys()
        free_energy = {k: 0.0 for k in lambda_names}
        var = {k: 0.0 for k in lambda_names}

        # iterate over lambda dimensions
        for i, lam_name in enumerate(lambda_names):
            lambda_array = getattr(lambdas, lam_name)

            dlam = get_lam_diff(lambda_array)

            # if dlam is not zero then bootstrap the replicas to get avg and var
            avg_var_by_state = np.array([compute_bs_error(state_data) if x != 0.0 else [0.0, 0.0] for
                                         state_data, x in zip(data[i], dlam)])

            avg_by_state = np.array([x[0] for x in avg_var_by_state])
            var_by_state = np.array([x[1] for x in avg_var_by_state])

            free_energy[lam_name] = np.trapz(avg_by_state, lambda_array)
            var[lam_name] = np.sum(np.square(dlam) * var_by_state)

        return free_energy, var

    def plot_du_by_dl(self):
        '''
        Function to plot dU/dlam vs state with include calculations for ci

        :return: None
        '''
        print('Plotting du/dlam graph...')
        avg_over_iterations = np.average(self.data, axis=3).transpose()
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.25 * 2, 4.5 * 2))

        keys = self.lambdas.schedule[0].keys()
        for state_data, lam_key in zip(avg_over_iterations, keys):
            lambda_array = getattr(self.lambdas, lam_key)
            # if lambdas not all zero or 1 then plot line
            if not all(v == 0 for v in lambda_array) and not all(v == 1 for v in lambda_array):
                # iteration over state
                avg_var_by_state = np.array([compute_bs_error(data) for data in state_data])
                avg_by_state = np.array([x[0] for x in avg_var_by_state])
                var_by_state = np.array([x[1] for x in avg_var_by_state])
                std_by_state = np.sqrt(var_by_state)

                dlam = get_lam_diff(lambda_array)
                x_vals = []
                y_vals = []
                err_vals = []
                for x, (y, std, dl) in enumerate(zip(avg_by_state, std_by_state, dlam)):
                    if dl != 0:
                        x_vals.append(x + 1)
                        y_vals.append(y)
                        err_vals.append(std)
                x_vals = np.array(x_vals)
                y_vals = np.array(y_vals)
                err_vals = np.array(err_vals)

                # avg_by_state = np.array([x if y != 0.0 else 0.0 for x, y in zip(avg_by_state, dlam)])
                # std_by_state = np.array([x if y != 0.0 else 0.0 for x, y in zip(std_by_state, dlam)])
                ax.set_xticks(range(1, self.shape[1]+1))
                ax.plot(x_vals, y_vals, label=lam_key, zorder=20)
                ax.fill_between(x_vals, y_vals + err_vals, y_vals - err_vals, alpha=0.15, zorder=1)
                ax.set_xlabel('State')
                ax.set_ylabel('$\delta U/ \delta \lambda$ [kcal/mol]')

        ax.set_xlim(1, self.shape[1])
        ax.legend(loc='upper left', bbox_to_anchor=(1.0, 1.025), labelspacing=1.2)
        plt.tight_layout()
        plt.savefig(os.path.join(self.analysis_dir, 'dUdlam'))
        plt.close('all')

def get_lam_diff(lambda_array):
    '''
    Compute difference between adjacent lambdas

    :param lambda_array: list of ints
    :return: numpy array for differences in input
    '''
    lambda_mp = np.array(np.multiply(0.5, np.add(lambda_array[1:], lambda_array[:-1])))
    lambda_mp = np.insert(lambda_mp, 0, lambda_array[0])
    lambda_mp = np.append(lambda_mp, lambda_array[-1])

    dlam = np.array(np.subtract(lambda_mp[1:], lambda_mp[:-1]))
    return dlam


def compute_bs_error(replicas):
    '''
    compute bootstrapped average and variance

    :param replicas: list, values of average dU/dlam in each replica

    :return: turple of floats, average and var of boot strapped result
    '''
    if len(replicas) > 1:
        boot = [resample(replicas, replace=True, n_samples=len(replicas)) for x in range(5000)]
        avg_per_boot_sample = np.average(boot, axis=1)
    else:
        return replicas[0], 0.0

    return np.average(avg_per_boot_sample), np.var(avg_per_boot_sample, ddof=1)




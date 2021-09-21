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
import math

class Lambdas(object):
    '''
    Class containing lambda schedule

    :param elec_edges: list, values for what fraction of the lambda schedule the electric potential should start and stop appearing
    :param ster_edges: list, values for what fraction of the lambda schedule the steric potential should start and stop appearing
    :param global_lambs: list of floats for global lambda schedule
    :param custom: dict, containing keys lambda_(sterics/electrostatics)_(appear/disappear) with list of floats for values
    :param debug: boolean, whether to print the schedule to terminal
    '''
    def __init__(self, elec_edges, ster_edges, global_lambs, custom=None, debug=True):

        if custom is None:

            self.lambda_sterics_appear = get_line(ster_edges[0], ster_edges[1], global_lambs, appear=True)
            self.lambda_electrostatics_appear = get_line(elec_edges[0], elec_edges[1], global_lambs, appear=True)

            self.lambda_sterics_disappear = get_line(1-ster_edges[1], 1-ster_edges[0], global_lambs, appear=False)
            self.lambda_electrostatics_disappear = get_line(1-elec_edges[1], 1-elec_edges[0], global_lambs, appear=False)

        else:
            self.lambda_sterics_appear = custom['lambda_sterics_appear']
            self.lambda_electrostatics_appear = custom['lambda_electrostatics_appear']
            self.lambda_sterics_disappear = custom['lambda_sterics_disappear']
            self.lambda_electrostatics_disappear = custom['lambda_electrostatics_disappear']

        if debug:
            print('Lambda sterics disappear: {}'.format(list(self.lambda_sterics_disappear)))
            print('Lambda sterics appear: {}'.format(list(self.lambda_sterics_appear)))
            print('Lambda electrostatics disappear: {}'.format(list(self.lambda_electrostatics_disappear)))
            print('Lambda electrostatics appear: {}'.format(list(self.lambda_electrostatics_appear)))


        assert (len(self.lambda_sterics_appear) == len(self.lambda_electrostatics_appear))
        assert (len(self.lambda_sterics_disappear) == len(self.lambda_electrostatics_appear))
        assert (len(self.lambda_sterics_disappear) == len(self.lambda_electrostatics_disappear))

        self.schedule = []
        for i, (su, sd, eu, ed) in enumerate(zip(self.lambda_sterics_appear, self.lambda_sterics_disappear,
                                                 self.lambda_electrostatics_appear, self.lambda_electrostatics_disappear)):
                param_vals = {'lambda_sterics_appear': su, 'lambda_sterics_disappear': sd,
                              'lambda_electrostatics_appear': eu, 'lambda_electrostatics_disappear': ed}
                self.schedule.append(param_vals)

    def update_attrs_from_schedule(self):
        '''
        Function to update the lambda attributes e.g. self.lambda_sterics_appear from self.schedule if it has changed.
        '''
        self.lambda_sterics_appear = [x['lambda_sterics_appear'] for x in self.schedule]
        self.lambda_electrostatics_appear = [x['lambda_electrostatics_appear'] for x in self.schedule]
        self.lambda_sterics_disappear = [x['lambda_sterics_disappear'] for x in self.schedule]
        self.lambda_electrostatics_disappear = [x['lambda_electrostatics_disappear'] for x in self.schedule]


def get_line(start, stop, global_lam, appear):
    '''
    Function to create a numpy array interpolating between end states to use in lambda schedule

    :param start: float, value at which line should begin to change
    :param stop:  float, value at which line should be finished changing
    :param global_lam: list, of values in the global lambdas
    :param appear: bool, True for appearing line i.e 0-> 1 False for disappear i.e 1->0

    :return: np.array([]), values for the line across windows

    '''

    if appear:
        line = [appear_func(start, stop, x) for x in global_lam]
    else:
        line = [disappear_func(start, stop, x) for x in global_lam]

    return np.array(line)


def appear_func(start, stop, x):
    '''
    Function to get the evaluation of a line constructed from two points (start, 0) and (stop, 1)

    :param start: float, when line should start appearing
    :param stop: float, when line should finish appearing
    :param x: float, value for which line should be evaluate at

    :return: float, y value of evaluated line

    '''
    if x < start:
        return 0
    if start <= x <= stop:
        dy = 1
        dx = stop - start
        return round((dy / dx) * (x - start) + 0, 6)
    else:
        return 1


def disappear_func(start, stop, x):
    '''
    Function to get the evaluation of a line constructed from two points (start, 1) and (stop, 0)

    :param start: float, when line should start disappearing
    :param stop: float, when line should finish disappearing
    :param x: float, value for which line should be evaluate at

    :return: float, y value of evaluated line

    '''
    if x < start:
        return 1
    if start <= x <= stop:
        dy = -1
        dx = stop - start
        return round((dy / dx) * (x - start) + 1, 6)
    else:
        return 0

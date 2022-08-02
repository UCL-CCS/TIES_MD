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

class Config():
    '''
    Class that reads the config file and set options.

    :param analysis_cfg: path to file containing general analysis params
    '''
    def __init__(self, analysis_cfg='./analysis.cfg'):
        general_args = read_config(analysis_cfg)
        self.engines = general_args['engines']

        #check all the config file args we need are present
        args_list = ['output_dir', 'temperature', 'legs', 'engines', 'data_root', 'exp_data',
                     'windows_mask', 'methods', 'distributions', 'rep_convg', 'sampling_convg',
                     'vdw_a', 'ele_a', 'vdw_d', 'ele_d']

        optional_args = ['namd_version']

        # check we have all required arguments
        for argument in args_list:
            if argument not in general_args.keys():
                raise ValueError('Missing option {} in configuration file'.format(argument))

        # check we have no unexpected arguments
        for argument in general_args.keys():
            if argument not in args_list+optional_args:
                raise ValueError('Argument {} not supported for this engine or at all.'
                                 ' Please remove from the TIES.cfg.'.format(argument))

        self.all_args = args_list+optional_args

        self.temperature = float(general_args['temperature'][0]) #unit = kelvin
        self.distributions = bool(int(general_args['distributions'][0]))
        self.engines_to_init = [x.lower() for x in general_args['engines']]
        self.simulation_legs = general_args['legs']
        self.data_root = general_args['data_root'][0]
        self.output_dir = general_args['output_dir'][0]
        self.namd_version = general_args['namd_version'][0]
        self.methods = general_args['methods']
        self.legs = general_args['legs']
        self.vdw_a = [float(x) for x in general_args['vdw_a']]
        self.vdw_d = [float(x) for x in general_args['vdw_d']]
        self.ele_a = [float(x) for x in general_args['ele_a']]
        self.ele_d = [float(x) for x in general_args['ele_d']]

        if self.simulation_legs == ['EDITME']:
            raise ValueError('Please set legs option in analysis.cfg')

        if general_args['windows_mask'][0] != 'None':
            self.windows_mask = [int(x) for x in general_args['windows_mask']]
        else:
            self.windows_mask = None

        if general_args['rep_convg'][0] != 'None':
            self.rep_convg = [int(x) for x in general_args['rep_convg']]
        else:
            self.rep_convg = None

        if general_args['sampling_convg'][0] != 'None':
            self.sampling_convg = [int(x) for x in general_args['sampling_convg']]
        else:
            self.sampling_convg = None

        if len(self.engines) == 0:
            raise ValueError('No support engines requested. Please choose from (NAMD2/NAMD3/OpenMM/GROMACS)')

        #process exp data
        # EXP data specifies what proteins and ligands to look at
        self.exp_data = eval(open(general_args['exp_data'][0]).read())

    def get_options(self):
        '''
        Function to print options stored in config class

        :return: None
        '''
        for arg in self.all_args:
            print('{}: {}'.format(arg, self.__getattribute__(arg)))


def read_config(config_file):
    '''
    Function that is reading config file and build dict of options.

    :param config_file: str, path to a config file

    :return: dict, containing all options to run analysis
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
                data = [s.replace("\t", "") for s in data]
                args_dict[data[0]] = data[1]
    args_dict = {k: v.split(',') for k, v in args_dict.items()}

    return args_dict
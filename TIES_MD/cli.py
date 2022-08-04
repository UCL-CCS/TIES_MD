#!/usr/local/bin/env python

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

from .TIES import TIES

from docopt import docopt
import os

# =============================================================================================
# COMMAND-LINE INTERFACE
# =============================================================================================

usage = """
TIES_MD
Command line input should be used as follows...
Usage:
TIES_MD [--devices=LIST] [--run_type=STRING] [--config_file=STRING] [--replica_mask=INT] [--windows_mask=INT] [--exp_name=STR]...
"""

def main(argv=None):
    '''
    Entry point for the TIES program.

    :param argv: dict, containing command line arguments
    '''
    args = docopt(usage, argv=argv, options_first=True)

    msg = 'No {0} specified using default {1}'

    #deal with command line options common for all engiens
    if args['--config_file']:
        config_file = args['--config_file']
        input_folder = config_file.split('/')[:-1]
        input_folder = '/'.join(input_folder)
    else:
        input_folder = os.getcwd()
        config_file = input_folder+'/TIES.cfg'
        print(msg.format('configuration file', config_file))

    run_types = ['run', 'setup', 'class']
    if args['--run_type']:
        run_type = args['--run_type']
        if run_type not in run_types:
            raise ValueError('Unsupported run type {}. Please choose from {}'.format(run_type, run_types))
    else:
        run_type = 'run'
        print(msg.format('Run type', run_type))

    if args['--exp_name']:
        exp_name = args['--exp_name'][0]
    else:
        exp_name = 'complex'
        print(msg.format('Experiment name', exp_name))

    if args['--replica_mask']:
        replica_mask = args['--replica_mask']
        replica_mask = replica_mask.split(',')
        replica_mask = [int(x) for x in replica_mask]
    else:
        replica_mask = None

    if args['--windows_mask']:
        window_masks = args['--windows_mask']
        window_masks = window_masks.split(',')
        window_masks = [int(x) for x in window_masks]
    else:
        window_masks = None

    # Read config file
    args_dict = read_config(config_file)

    #openmm specific command line options
    not_openmm_msg = 'Option {} only supported for OpenMM TIES' 
    if args_dict['engine'].lower() != 'openmm':
        not_openmm = True
    else:
        not_openmm = False

    if args['--devices']:
        if not_openmm:
            raise ValueError(not_openmm_msg.format('--devices'))
        devices = args['--devices']
        devices = devices.split(',')
        devices = [int(x) for x in devices]
    else:
        devices = None

    #removed this as an option there is no need to expose it for now
    periodic = True

    TIES(input_folder, exp_name, run_type, devices, replica_mask, window_masks, **args_dict)

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

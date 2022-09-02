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

from .TIES import TIES, read_config

from docopt import docopt
import os

# =============================================================================================
# COMMAND-LINE INTERFACE
# =============================================================================================

usage = """
TIES_MD
Command line input should be used as follows...
Usage:
TIES_MD [--devices=LIST] [--run_type=STRING] [--config_file=STRING] [--rep_id=INT] [--windows_mask=LIST] [--exp_name=STR]...
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
        input_folder = os.path.abspath('/'.join(input_folder))
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

    if args['--windows_mask']:
        mask = args['--windows_mask']
        mask = mask.split(',')
        mask = [int(x) for x in mask]
    else:
        mask = None

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

    if args['--rep_id']:
        if not_openmm:
            raise ValueError(not_openmm_msg.format('--rep_id'))
        rep_id = args['--rep_id']
        rep_id = int(rep_id)
    else:
        rep_id = None
        print(msg.format('node id string', 'None'))

    # removed this as an option there is no need to expose it for now
    periodic = True

    TIES(input_folder, exp_name, run_type, devices, rep_id, mask, periodic, **args_dict)


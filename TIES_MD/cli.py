#!/usr/local/bin/env python

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
TIES_MD [--devices=LIST] [--run_type=STRING] [--config_file=STRING] [--periodic=BOOL] [--node_id=STRING] [--windows_mask=LIST] [--exp_name=STR]...
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
        devices = [0]
        print(msg.format('devices', devices))

    if args['--node_id']:
        if not_openmm:
            raise ValueError(not_openmm_msg.format('--node_id'))
        node_id = args['--node_id']
        node_id = str(node_id)
    else:
        node_id = None
        #node_id will be set to _alpha in TIES.py after options check
        print(msg.format('node id string', '_alpha'))

    if args['--windows_mask']:
        if not_openmm:
            raise ValueError(not_openmm_msg.format('--windows_mask'))
        mask = args['--windows_mask']
        mask = mask.split(',')
        mask = [int(x) for x in mask]
    else:
        mask = None

    if args['--periodic']:
        if not_openmm:
            raise ValueError(not_openmm_msg.format('--periodic'))
        periodic = bool(int(args['--periodic']))
    else:
        periodic = True
        print(msg.format('spatial periodicity', periodic))

    TIES(input_folder, run_type, exp_name, devices, node_id, mask, periodic, args_dict)

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

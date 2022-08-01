#!/usr/bin/env python

import os
import argparse
import glob

from .config import Config
import numpy as np
from .engines.openmm import OpenMM
from .engines.namd import NAMD
from .engines._numpy import Numpy

class Analysis():
    '''
    Class that is setting up all engines and methods for analysis.

    :param cfg: dict, holding all options
    '''
    def __init__(self, cfg):
        self.cfg = cfg
        # Process engines
        self.engines = []
        if 'namd2' in cfg.engines_to_init or 'namd3' in cfg.engines_to_init:
            for method in cfg.methods:
                namd = NAMD(method.upper(), cfg.output_dir, cfg.windows_mask, cfg.distributions, cfg.rep_convg,
                            cfg.sampling_convg, cfg.vdw_a, cfg.vdw_d, cfg.ele_a, cfg.ele_d, cfg.namd_version)
                self.engines.append(namd)

        if 'openmm' in cfg.engines_to_init:
            for method in cfg.methods:
                openmm = OpenMM(method.upper(), cfg.output_dir, cfg.windows_mask, cfg.distributions, cfg.rep_convg,
                                cfg.sampling_convg, cfg.vdw_a, cfg.vdw_d, cfg.ele_a, cfg.ele_d, fep_combine_reps=False)
                self.engines.append(openmm)

        if 'numpy' in cfg.engines_to_init:
            for method in cfg.methods:
                numpy = Numpy(method.upper(), cfg.output_dir, cfg.windows_mask, cfg.distributions, cfg.rep_convg,
                              cfg.sampling_convg, cfg.vdw_a, cfg.vdw_d, cfg.ele_a, cfg.ele_d, fep_combine_reps=False)
                self.engines.append(numpy)

        if len(self.engines) == 0:
            raise ValueError('No valid engine found. Select from namd2, namd3 or openmm')

    def run(self):
        '''
        Function to run the analysis.

        :return: None
        '''
        result = {}
        for engine in self.engines:
            engine_id = engine.name+'_'+engine.method
            if engine_id not in result:
                result[engine_id] = {}
            else:
                raise ValueError('Repeated engine name {}'.format(engine_id))

            for prot_name, prot_data in self.cfg.exp_data.items():
                if prot_name not in result[engine_id]:
                    result[engine_id][prot_name] = {}
                else:
                    raise ValueError('Repeated protein name {}'.format(prot_name))

                for lig_name, lig_data in prot_data.items():
                    if lig_name not in result[engine_id][prot_name]:
                        result[engine_id][prot_name][lig_name] = {}
                    else:
                        raise ValueError('Repeated ligand name {}'.format(lig_name))

                    nice_print('{} {} {}'.format(engine_id, prot_name, lig_name))
                    leg_results = {leg: 0.00 for leg in self.cfg.simulation_legs}

                    for leg in self.cfg.simulation_legs:
                        nice_print(leg)
                        leg_results[leg] = engine.run_analysis(self.cfg.data_root, self.cfg.temperature,
                                                               prot_name, lig_name, leg)

                    if len(self.cfg.simulation_legs) == 2:
                        print('Two thermodynamic legs found assuming this is a ddG calculation')
                        leg1 = self.cfg.simulation_legs[0]
                        leg2 = self.cfg.simulation_legs[1]
                        print('')
                        print('{} dG = {}: SEM = {}'.format(leg1, leg_results[leg1][0], leg_results[leg1][1]))
                        print('{} dG = {}: SEM = {}'.format(leg2, leg_results[leg2][0], leg_results[leg2][1]))
                        print('Computing {} - {}'.format(leg1, leg2))

                        ddg = leg_results[leg1][0] - leg_results[leg2][0]
                        ddg_err = np.sqrt(np.square(leg_results[leg1][1]) + np.square(leg_results[leg2][1]))
                        print('')
                        print('ddG = {}: SEM = {}'.format(ddg, ddg_err))
                        result[engine_id][prot_name][lig_name] = [ddg, ddg_err]

                    elif len(self.cfg.simulation_legs) == 1:
                        leg1 = self.cfg.simulation_legs[0]
                        dg = leg_results[leg1][0]
                        dg_err = leg_results[leg1][1]
                        print('')
                        print('dG = {}: SEM = {}'.format(dg, dg_err))
                        result[engine_id][prot_name][lig_name] = [dg, dg_err]

                    else:
                        for leg in self.cfg.simulation_legs:
                            print('result = {} SEM = {} for leg {}'.format(*leg_results[leg], leg))
                        fin_result = {leg: leg_results[leg] for leg in self.cfg.simulation_legs}
                        result[engine_id][prot_name][lig_name] = fin_result

        with open('./result.dat', 'w') as f:
            print(result, file=f)


def nice_print(string):
    '''
    Function to print to terminal with padding.

    :param string: str to format

    :return: None
    '''
    string = string.center(75, '#')
    print(string)


def main():
    '''
    Main function where we chose to setup or run analysis.

    :return: None
    '''
    nice_print('TIES Analysis')

    parser = argparse.ArgumentParser(description='Program for the analysis of TIES_MD outputs.')
    parser.add_argument("--run_type", type=str, help="What actions to perform [setup/run]", required=False)
    args = parser.parse_args()

    if args.run_type == 'run' or args.run_type is None:
        cfg = Config()
        print('Running with options:')
        cfg.get_options()
        ana = Analysis(cfg)
        ana.run()
        nice_print('END')
    elif args.run_type == 'setup':
        print('Generating dummy exp data...')
        print('If legs have been set in analysis.cfg you can now run TIES_ana to get results.')
        #looking for results files
        populated_openmm_results = list(glob.iglob(os.path.join('*', '*', '*', 'LAMBDA*', 'rep*', 'results', '*.npy')))
        populated_namd_results = list(glob.iglob(os.path.join('*', '*', '*', 'LAMBDA*', 'rep*', 'simulation', '*.alch')))
        populated_openmm_results = list(set([x.split('LAMBDA')[0] for x in populated_openmm_results]))
        populated_namd_results = list(set([x.split('LAMBDA')[0] for x in populated_namd_results]))
        if len(populated_namd_results) > 0 and len(populated_openmm_results) > 0:
            raise ValueError('Mixed OpenMM and NAMD results found please conduct analysis separately for each engine')
        pop_results = populated_openmm_results+populated_namd_results
        exp_dat = {}
        for x in pop_results:
            _data = x.split(os.sep)
            ligand = {_data[1]: [0.00, 0.00]}

            if _data[0] in exp_dat:
                if _data[1] in exp_dat[_data[0]]:
                    #ignore if we have this included already
                    continue
                else:
                    #if we have seen this system befor just add the new ligand
                    exp_dat[_data[0]].update(ligand)
            else:
                #if we have never see this in any form add all info
                exp_dat.update({_data[0]: ligand})

        with open('exp.dat', 'w') as f:
            print(exp_dat, file=f)

        nice_print('END')
    else:
        raise ValueError('Unknown value for arg run_type. Select from [setup/run]')




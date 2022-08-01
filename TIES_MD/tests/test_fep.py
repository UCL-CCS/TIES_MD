#!/usr/bin/env python

import unittest
import numpy as np

from TIES_MD.ties_analysis.methods.FEP import *
from TIES_MD.ties_analysis.engines.openmm import Lambdas

from pymbar.testsystems.harmonic_oscillators import HarmonicOscillatorsTestCase
from pymbar import MBAR, timeseries

#GLOBAL CONSTANTS
KBT = 0.5961619840741075 #unit kilocalorie/mole assumed 300k


class Test_MBAR(unittest.TestCase):

        def test_mbar(self):

            harm = HarmonicOscillatorsTestCase()
            x_kn, u_kln, N_k = harm.sample(N_k=[1000, 1000, 1000, 1000, 1000], mode='u_kln')

            #These lambda values are dummy just to init the class MBAR_Analysis
            lam = Lambdas(vdw_a=[0.0, 0.5, 1.0, 1.0, 1.0],
                          vdw_d=[0.0, 0.0, 0.3333333333333333, 0.6666666666666667, 1.0],
                          ele_a=[1.0, 1.0, 0.6666666666666667, 0.33333333333333337, 0.0],
                          ele_d=[1.0, 0.5, 0.0, 0.0, 0.0])

            #use TIES code
            mbar = MBAR_Analysis(MBARs=np.array([u_kln]), temp=300, lambdas=lam, analysis_dir=None, decorrelate=True)
            avg, var = mbar.analysis()

            #check against minimal MBAR example
            nstates = 5
            N_k = np.zeros([nstates], np.int32)  # number of uncorrelated samples
            for k in range(nstates):
                [_, g, __] = timeseries.detect_equilibration(u_kln[k, k, :])
                indices = timeseries.subsample_correlated_data(u_kln[k, k, :], g=g)
                N_k[k] = len(indices)
                u_kln[k, :, 0:N_k[k]] = u_kln[k, :, indices].T
            mbar = MBAR(u_kln, N_k)
            mbar_results = mbar.compute_free_energy_differences(return_theta=True)

            result = ([mbar_results['Delta_f'][0, len(N_k) - 1] * KBT,
                       mbar_results['dDelta_f'][0, len(N_k) - 1] * KBT])  # units kilocalorie per mol

            self.assertEqual(round(avg, 5), round(result[0], 5), 'MBAR analysis failed')

            #this is the check carried out in chodera code to test against analytic resutls.
            ans = harm.analytical_free_energies()[-1]*KBT
            z = (result[0] - ans) / result[1]
            self.assertEqual(round(z/12, 0), 0, 'MBAR analysis failed')


if __name__ == '__main__':
    unittest.main()
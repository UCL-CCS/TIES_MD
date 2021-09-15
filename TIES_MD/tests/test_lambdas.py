#!/usr/bin/env python

import unittest
from TIES_MD.lambdas import *

class Test_Lambdas(unittest.TestCase):

    def test_lambda(self):
        #test we can construct the lambda schedual from 2017 paper: https://doi.org/10.1021/acs.jctc.6b00979
        test_vals = [[0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0],
                     [1.0, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.0],
                     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.090909, 0.272727, 0.454545, 0.636364, 0.818182, 0.909091, 1.0],
                     [1.0, 0.909091, 0.818182, 0.636364, 0.454545, 0.272727, 0.090909, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]

        lam = Lambdas([0.45, 1.0], [0.0, 1.0], [0.00, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.00],
                      debug=False)

        sa = [x['lambda_sterics_appear'] for x in lam.schedule]
        sd = [x['lambda_sterics_disappear'] for x in lam.schedule]
        ea = [x['lambda_electrostatics_appear'] for x in lam.schedule]
        ed = [x['lambda_electrostatics_disappear'] for x in lam.schedule]

        print(sa)
        print(sd)
        print(ea)
        print(ed)

        self.assertEqual(sa, test_vals[0], 'lambda_sterics_appear incorrectly constructed')
        self.assertEqual(sd, test_vals[1], 'lambda_sterics_disappear incorrectly constructed')
        self.assertEqual(ea, test_vals[2], 'lambda_electrostatics_appear incorrectly constructed')
        self.assertEqual(ed, test_vals[3], 'lambda_electrostatics_disappear incorrectly constructed')

        del lam.schedule[0]
        del lam.schedule[10]

        del test_vals[0][0]
        del test_vals[1][0]
        del test_vals[2][0]
        del test_vals[3][0]

        del test_vals[0][10]
        del test_vals[1][10]
        del test_vals[2][10]
        del test_vals[3][10]

        lam.update_attrs_from_schedule()

        sa = [x['lambda_sterics_appear'] for x in lam.schedule]
        sd = [x['lambda_sterics_disappear'] for x in lam.schedule]
        ea = [x['lambda_electrostatics_appear'] for x in lam.schedule]
        ed = [x['lambda_electrostatics_disappear'] for x in lam.schedule]

        self.assertEqual(sa, test_vals[0], 'lambda_sterics_appear incorrectly updated')
        self.assertEqual(sd, test_vals[1], 'lambda_sterics_disappear incorrectly updated')
        self.assertEqual(ea, test_vals[2], 'lambda_electrostatics_appear incorrectly updated')
        self.assertEqual(ed, test_vals[3], 'lambda_electrostatics_disappear incorrectly updated')


if __name__ == '__main__':
    unittest.main()


#!/usr/bin/env python

import unittest
import numpy as np

from TIES_MD.ties_analysis.methods.TI import *
from TIES_MD.ties_analysis.engines.openmm import Lambdas

class Test_TI(unittest.TestCase):

    def test_bootstrapping(self):

        mus = [1, 3, -0.3]
        sigmas  = [0.02, 1.2, 1.5]


        for mu, s in zip(mus, sigmas):

            data = np.random.normal(mu, s, 5000)

            sem = s/np.sqrt(len(data))
            mean = (np.mean(data))

            #computed with code
            avg, var = compute_bs_error(data)

            self.assertEqual(round(avg, 2), round(mean, 2), 'Mean calculation failed')
            self.assertEqual(round(np.sqrt(var), 2), round(sem, 2), 'SEM calculation failed')


    def test_intergration(self):

        def f(x):
            return x

        def f1(x):
            return -x

        def f2(x):
            return x*2 + 4

        def f3(x):
            return x**2 - 0.88

        #using of odd number of steps can confuse range generation and it may not run 0->10 as intended
        steps = [2, 2, 50, 500]
        answers = [200, -200, 560.0, 1298.132]

        for true_ans, function, s in zip(answers, [f, f1, f2, f3], steps):
            x,y = get_data(function, s)
            data = np.zeros([1, len(y), 4, 1])
            y = np.array([[[r] for r in y]])
            for i in range(4):
                data[:, :, i, :] = y

            lam = Lambdas(x, x, x, x)
            ti = TI_Analysis(data, lam, None)
            ans, err = ti.analysis()

            self.assertEqual(round(ans, 2), round(true_ans, 2), 'Integration failed for {}'.format(str(function)))


def get_data(f, step):
    result  = []
    for x in np.arange(0, 10+(1/step), 1/step):
        result.append(f(x))
    return np.arange(0, 10+(1/step), 1/step), result


if __name__ == '__main__':
    unittest.main()


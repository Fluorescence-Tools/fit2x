from __future__ import division
import unittest

import numpy as np
import fit2x

np.random.seed(42)

x = np.arange(0, 1024)
ym = np.ones_like(x) * 500
yd = np.random.poisson(ym)
data = fit2x.DecayCurve(x, yd)
model = fit2x.DecayCurve(x, ym)
score_settings = {
    "model": model,
    "data": data,
    "start": 0,
    "stop": -1
}


class Tests(unittest.TestCase):

    def test_score_1(self):
        dc = fit2x.DecayScore(**score_settings)
        self.assertAlmostEqual(dc.score, 1114.2359946704964)
        dc.model.y = ym * 0.95
        self.assertAlmostEqual(dc.score, 2185.4654851901237)

    def test_score_range(self):
        dc = fit2x.DecayScore(**score_settings)
        dc.start, dc.stop = 10, 20
        wr = dc.weighted_residuals
        ref = np.array([-0.13456839, -0.95951449,  0.26673253,  0.31088091, -0.58908623,
                        0.53033009, -1.38379681, -0.72727273, -0.58908623,  0.13376339])
        (np.allclose(wr, ref), True)

    def test_score_type(self):
        dc = fit2x.DecayScore(**score_settings)
        dc.model.y = ym

        dc.score_type = 'default'
        self.assertEqual(dc.score_type, 'default')
        self.assertAlmostEqual(dc.score, 1114.2359946704964)

        dc.score_type = 'poisson'
        self.assertEqual(dc.score_type, 'poisson')
        self.assertAlmostEqual(dc.score, 1103.3028342960324)

        dc.score_type = 'neyman'
        self.assertEqual(dc.score_type, 'neyman')
        self.assertAlmostEqual(dc.score, 1114.2359946704964)

        dc.score_type = 'pearson'
        self.assertEqual(dc.score_type, 'pearson')
        self.assertAlmostEqual(dc.score, 2.0540000000000003)

        dc.score_type = 'gauss'
        self.assertEqual(dc.score_type, 'gauss')
        self.assertAlmostEqual(dc.score, 1103.2583213870241)

        dc.score_type = 'cnp'
        self.assertEqual(dc.score_type, 'cnp')
        self.assertAlmostEqual(dc.score, 1104.466664890166)

    def test_setter_getter(self):
        dc = fit2x.DecayScore(**score_settings)
        self.assertEqual(np.allclose(dc.model.y, model.y), True)
        self.assertEqual(np.allclose(dc.data.y, data.y), True)

        x = np.arange(10, dtype=np.float64)
        a1 = np.ones_like(x, dtype=np.float64)
        a2 = np.ones_like(x, dtype=np.float64) * 2
        c1 = fit2x.DecayCurve(x, a1)
        c2 = fit2x.DecayCurve(x, a2)
        dc.data = c1
        dc.model = c2
        self.assertEqual(np.allclose(dc.data.y, c1.y), True)
        self.assertEqual(np.allclose(dc.model.y, c2.y), True)

    @unittest.expectedFailure
    def test_setter_getter_2(self):
        dc = fit2x.DecayScore(**score_settings)
        self.assertEqual(np.allclose(dc.model.y, model.y), True)
        self.assertEqual(np.allclose(dc.data.y, data.y), True)

        x = np.arange(10, dtype=np.float64)
        a1 = np.ones_like(x, dtype=np.float64)
        a2 = np.ones_like(x, dtype=np.float64) * 2
        # Works on pointers reference counting issues
        # see: https://stackoverflow.com/questions/20029377/why-does-swig-appear-to-corrupt-contents-of-a-member-class
        dc.data = fit2x.DecayCurve(x, a1)
        dc.model = fit2x.DecayCurve(x, a2)
        self.assertEqual(np.allclose(dc.data.y, a1), True)
        self.assertEqual(np.allclose(dc.model.y, a2), True)


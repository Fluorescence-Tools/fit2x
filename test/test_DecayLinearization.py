from __future__ import division
import unittest

import numpy as np
import fit2x


class Tests(unittest.TestCase):

    def test_lin(self):
        x = np.arange(0, 16)
        y = np.ones_like(x) * 10000
        model = fit2x.DecayCurve(x, y)

        lin = np.sin(x) * 0.2 + 0.8
        linearization = fit2x.DecayCurve(x, lin)
        dl = fit2x.DecayLinearization(
            linearization_table=linearization,
            start=0, stop=-1, active=True
        )
        dl.add(model)

        ref = np.array([8000., 9682.94196962, 9818.59485365, 8282.24001612,
                        6486.39500938, 6082.15145067, 7441.1690036, 9313.97319744,
                        9978.71649325, 8824.23697048, 6911.95777822, 6000.0195869,
                        6926.854164, 8840.33407365, 9981.21471139, 9300.57568031])
        self.assertEqual(np.allclose(model.y, ref), True)

        model = fit2x.DecayCurve(x, y)
        dl.active = False
        dl.add(model)
        self.assertEqual(np.allclose(model.y, y), True)


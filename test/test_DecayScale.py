from __future__ import division
import unittest

import numpy as np
import fit2x

x = np.arange(0, 16)
y = np.ones_like(x) * 100
y2 = y * 10
data = fit2x.DecayCurve(x, y2)
settings = {
    "data": data,
    "constant_background": 0,
    "active": True,
    "start": 0,
    "stop": -1
}


class Tests(unittest.TestCase):

    def test_scale_init(self):
        ds = fit2x.DecayScale(**settings)
        self.assertEqual(np.allclose(ds.data.y, data.y), True)
        self.assertEqual(np.allclose(ds.data.x, data.x), True)
        self.assertEqual(ds.constant_background, 0.0)
        self.assertEqual(ds.active, True)

    def test_scale_scale_1(self):
        x = np.arange(0, 16)
        y = np.ones_like(x) * 100
        y2 = y * 10
        data = fit2x.DecayCurve(x, y2)
        settings = {
            "data": data,
            "constant_background": 0,
            "active": True,
            "start": 0,
            "stop": -1
        }
        ds = fit2x.DecayScale(**settings)
        model = fit2x.DecayCurve(x, y)
        ds.add(model)

        self.assertEqual(np.allclose(model.y, data.y), True)

        ds.active = False
        self.assertEqual(ds.active, False)

        model = fit2x.DecayCurve(x, y)
        ds.add(model)
        self.assertEqual(np.allclose(model.y, y), True)

    def test_scale_scale_2(self):
        import numpy as np
        import fit2x

        x = np.arange(0, 16)
        y = np.ones_like(x) * 100
        y2 = y * 10
        data = fit2x.DecayCurve(x, y2)
        settings = {
            "data": data,
            "constant_background": 0,
            "active": True,
            "start": 0,
            "stop": -1
        }

        ds = fit2x.DecayScale(**settings)
        model = fit2x.DecayCurve(x, y)
        self.assertEqual(np.allclose(ds.data.y, np.ones_like(ds.data.y) * 1000), True)
        self.assertEqual(np.allclose(model.y, np.ones_like(ds.data.y) * 100), True)
        ds.add(model)
        self.assertEqual(np.sum(model.y), np.sum(ds.data.y))

    def test_scale_scale_3(self):
        ds = fit2x.DecayScale(**settings)
        # Scaling considering noise & background of data
        ds.constant_background = 100
        # Set number of photons to negative value for scaling
        model = fit2x.DecayCurve(x, y)
        ds.add(model)
        self.assertEqual(np.allclose(model.y, ds.data.y - 100), True)
        self.assertEqual(ds.number_of_photons, 16000)

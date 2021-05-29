from __future__ import division
import unittest

import numpy as np
import fit2x


class Tests(unittest.TestCase):

    def test_pileup(self):
        x = np.arange(0, 16)
        y = np.ones_like(x) * 10000
        decay_settings = {
            "x": x,
            "y": y,
            "acquisition_time": 0.1
        }
        data = fit2x.DecayCurve(**decay_settings)
        model = fit2x.DecayCurve(**decay_settings)

        pileup_settings = {
            "data": data,
            "pile_up_model": "coates",
            "repetition_rate": 100.0,
            "instrument_dead_time": 120.0,
            "active": True
        }
        dp = fit2x.DecayPileup(**pileup_settings)
        self.assertEqual(dp.repetition_rate, 100.0)
        self.assertEqual(dp.instrument_dead_time, 120.0)

        dp.add(model)
        ref = np.array([10093.8673586, 10081.35171081, 10068.83606301, 10056.32041521,
                        10043.8047674, 10031.28911959, 10018.77347177, 10006.25782395,
                        9993.74217613, 9981.2265283, 9968.71088047, 9956.19523264,
                        9943.6795848, 9931.16393695, 9918.64828911, 9906.13264125])
        self.assertEqual(np.allclose(model.y, ref), True)

    def test_pileup_setter_getter(self):
        x = np.arange(0, 16)
        y = np.ones_like(x) * 10000
        decay_settings = {
            "x": x,
            "y": y,
            "acquisition_time": 0.1
        }
        data = fit2x.DecayCurve(**decay_settings)
        model = fit2x.DecayCurve(**decay_settings)

        dp = fit2x.DecayPileup()

        dp.repetition_rate = 13.0
        self.assertEqual(dp.repetition_rate, 13.0)

        dp.pile_up_model = "coates"
        self.assertEqual(dp.pile_up_model, "coates")

        dp.instrument_dead_time = 89
        self.assertEqual(dp.instrument_dead_time, 89)

        dp.data = data
        (np.allclose(dp.data.x, data.x), True)
        (np.allclose(dp.data.y, data.y), True)

        dp.resize(10)
        self.assertEqual(len(dp.data.y), 10)

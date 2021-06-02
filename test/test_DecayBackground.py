from __future__ import division
import unittest

import numpy as np
import fit2x


class Tests(unittest.TestCase):

    def test_constant_background(self):
        x = np.arange(0, 4, dtype=np.float64)
        y = np.zeros_like(x)
        dc = fit2x.DecayCurve(x, y)
        self.assertEqual(np.allclose(dc.y, y), True)
        bg = fit2x.DecayPattern(5.0)
        bg.add(dc)
        self.assertEqual(np.allclose(dc.y, y + 5), True)
        bg.add(dc)
        self.assertEqual(np.allclose(dc.y, y + 10), True)
        bg.constant_offset = -100
        self.assertEqual(bg.constant_offset, -100)
        bg.add(dc)
        self.assertEqual(np.allclose(dc.y, y - 90), True)

    def test_pattern(self):
        x = np.linspace(0, 10, 8, dtype=np.float64)
        y = np.sin(x) + 1.0
        dc = fit2x.DecayCurve(x, y)
        (np.allclose(dc.y, y), True)

        bg_pattern = x
        bg_pattern_curve = fit2x.DecayCurve(x, bg_pattern)
        f = 0.8
        offset = 1.0
        bg = fit2x.DecayPattern(
            constant_offset=offset,
            pattern=bg_pattern_curve,
            pattern_fraction=f
        )
        bg.add(dc)
        r = y * (1-f) + bg_pattern * f * y.sum() / sum(bg_pattern) + offset
        self.assertEqual(np.allclose(r, dc.y), True)

from __future__ import division

import unittest
import numpy as np
import scipy.stats

import fit2x
from compute_irf import model_irf


irf, time_axis = model_irf(
    n_channels=64,
    period=32.,
    irf_position_p=2.0,
    irf_position_s=18.0,
    irf_width=0.25
)
bg = np.zeros_like(irf) + 0.2


class Tests(unittest.TestCase):

    def test_fit25(self):
        import fit2x
        import scipy.stats
        n_channels = 32
        irf_position_p = 2.0
        irf_position_s = 18.0
        irf_width = 0.25
        irf, time_axis = model_irf(
            n_channels=n_channels,
            period=32.,
            irf_position_p=irf_position_p,
            irf_position_s=irf_position_s,
            irf_width=irf_width
        )
        data = np.array(
            [0, 0, 0, 1, 9, 7, 5, 5, 5, 2, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 2, 2, 2, 3, 0, 1, 0,
            1, 1, 1, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        )
        bg = np.zeros_like(irf)
        settings = {
            'dt': time_axis[1] - time_axis[0],
            'g_factor': 1.0,
            'l1': 0.1,
            'l2': 0.1,
            'period': 32.0,
            'convolution_stop': 31,
            'irf': irf,
            'background': bg,
            'verbose': True
        }
        fit25 = fit2x.Fit25(**settings)
        # import pylab as p
        # p.plot(s)
        # p.plot(model_ref)
        # p.show()
        r0 = 0.38
        tau1, tau2, tau3, tau4 = 0.5, 1.0, 2.0, 4.0
        gamma = 0.02
        x = np.array([tau1, tau2, tau3, tau4, gamma, r0])
        fixed = np.array([0, 0, 0, 0, 1, 1])
        r = fit25(
            data=data,
            initial_values=x,
            fixed=fixed
        )
        best_tau = r['x'][0]
        self.assertEqual(best_tau, 2.0)

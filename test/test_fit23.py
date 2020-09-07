from __future__ import division

import unittest
import numpy as np

import fit2x
from compute_irf import model_irf


class Tests(unittest.TestCase):

    def test_modelf23(self):
        irf, time_axis = model_irf(
            n_channels=64,
            period=32.,
            irf_position_p=2.0,
            irf_position_s=18.0,
            irf_width=0.25
        )
        bg = np.zeros_like(irf) + 0.2
        tau, gamma, r0, rho = 2.0, 0.01, 0.38, 1.2
        np.random.seed(0)
        param = np.array([tau, gamma, r0, rho])
        period, g, l1, l2, conv_stop = 32, 1.0, 0.1, 0.1, 63
        conv_stop = min(len(time_axis), conv_stop)
        corrections = np.array([period, g, l1, l2, conv_stop])
        dt = time_axis[1] - time_axis[0]
        out = np.zeros_like(bg)
        fit2x.modelf23(param, irf, bg, dt, corrections, out)
        out_ref = np.array(
            [0.00201599, 0.00201501, 0.0020141, 0.00201328, 0.00201637,
             0.00214542, 0.00374779, 0.01133864, 0.02594058, 0.03749508,
             0.03937617, 0.03621302, 0.03262197, 0.02954715, 0.02694633,
             0.02472382, 0.02280411, 0.02112878, 0.01965254, 0.0183401,
             0.01716387, 0.01610213, 0.0151377, 0.01425687, 0.01344863,
             0.01270406, 0.01201585, 0.01137795, 0.01078531, 0.01023366,
             0.00971935, 0.00923921, 0.00879051, 0.00837082, 0.00797797,
             0.00761005, 0.00726529, 0.00694212, 0.0066391, 0.00635489,
             0.00608826, 0.0058381, 0.00560335, 0.00538303, 0.00517624,
             0.00498214, 0.00479994, 0.00462889, 0.00446831, 0.00431756,
             0.00417603, 0.00404314, 0.00391838, 0.00380125, 0.00369126,
             0.003588, 0.00349105, 0.00340001, 0.00331454, 0.00323428,
             0.00315893, 0.00308817, 0.00302174, 0.00295936, 0.00201549,
             0.00201455, 0.00201366, 0.00201303, 0.00202367, 0.00225593,
             0.00406495, 0.00976845, 0.01771813, 0.02272655, 0.02380904,
             0.02338395, 0.02265528, 0.02182821, 0.02094789, 0.02004244,
             0.01913224, 0.01823191, 0.01735177, 0.01649889, 0.01567793,
             0.01489176, 0.01414192, 0.01342897, 0.01275278, 0.01211273,
             0.01150783, 0.01093689, 0.01039852, 0.0098913, 0.00941373,
             0.00896432, 0.00854158, 0.00814407, 0.00777038, 0.00741918,
             0.00708916, 0.00677909, 0.0064878, 0.00621418, 0.00595718,
             0.00571581, 0.00548912, 0.00527623, 0.00507631, 0.00488858,
             0.00471228, 0.00454674, 0.0043913, 0.00424533, 0.00410828,
             0.00397958, 0.00385874, 0.00374527, 0.00363873, 0.00353869,
             0.00344476, 0.00335656, 0.00327375, 0.00319599, 0.00312297,
             0.00305442, 0.00299005, 0.00292961]
        )
        self.assertEqual(
            np.allclose(out, out_ref),
            True
        )
        # p.plot(time_axis, out)
        # p.show()

    def test_modelf23(self):
        # parameters used in different tests
        n_channels = 32
        n_photons = 60
        irf_position_p = 2.0
        irf_position_s = 18.0
        irf_width = 0.25
        period, g, l1, l2, conv_stop = 32, 1.0, 0.1, 0.1, 255
        tau, gamma, r0, rho = 2.0, 0.01, 0.38, 1.2

        # setup some parameters
        irf_np, time_axis = model_irf(
            n_channels=n_channels,
            period=period,
            irf_position_p=irf_position_p,
            irf_position_s=irf_position_s,
            irf_width=irf_width
        )
        dt = time_axis[1] - time_axis[0]
        conv_stop = len(time_axis) / 2 - 1
        param = np.array([tau, gamma, r0, rho])
        corrections = np.array([period, g, l1, l2, conv_stop])
        # compute a model function that is later used as "data"
        model = np.zeros_like(time_axis)
        bg = np.zeros_like(time_axis)
        fit2x.modelf23(param, irf_np, bg, dt, corrections, model)
        model_ref = np.array(
            [2.19319033e-08, 1.84413639e-08, 3.34901093e-05, 1.26870040e-02,
             9.39283373e-02, 1.22572899e-01, 9.04619781e-02, 6.48528685e-02,
             4.76001206e-02, 3.55488448e-02, 2.68771272e-02, 2.04940827e-02,
             1.57172794e-02, 1.21005483e-02, 9.34005838e-03, 7.22159597e-03,
             5.58990082e-03, 4.33007676e-03, 3.35581222e-03, 2.60158411e-03,
             2.01729163e-03, 1.56443990e-03, 1.21335530e-03, 9.41114763e-04,
             7.29984906e-04, 5.66234218e-04, 4.39223436e-04, 3.40705850e-04,
             2.64287576e-04, 2.05010412e-04, 1.59029028e-04, 1.23360965e-04,
             1.86812539e-08, 1.77852830e-07, 4.86291035e-04, 2.41480990e-02,
             5.93758154e-02, 6.60604491e-02, 5.44150218e-02, 4.37331561e-02,
             3.46979744e-02, 2.73087097e-02, 2.13834235e-02, 1.66888768e-02,
             1.29973503e-02, 1.01084281e-02, 7.85456846e-03, 6.09967166e-03,
             4.73504656e-03, 3.67479642e-03, 2.85148543e-03, 2.21239392e-03,
             1.71641883e-03, 1.33157038e-03, 1.03297999e-03, 8.01329509e-04,
             6.21619684e-04, 4.82208333e-04, 3.74060853e-04, 2.90167144e-04,
             2.25088434e-04, 1.74605311e-04, 1.35444470e-04, 1.05066633e-04]
        )
        # import pylab as p
        # p.plot(model)
        # p.show()
        print(model)
        self.assertEqual(
            np.allclose(model, model_ref), True
        )
        # add poisson noise to model and use as data
        np.random.seed(0)
        data = np.random.poisson(model * n_photons)
        data_ref = np.array(
            [0, 0, 0, 1, 9, 7, 5, 5, 5, 2, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 2, 2, 2, 3, 0, 1, 0,
             1, 1, 1, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        )
        self.assertEqual(
            np.allclose(data, data_ref), True
        )

    def test_target23(self):
        n_channels = 32
        irf_position_p = 2.0
        irf_position_s = 18.0
        irf_width = 0.25
        period, g, l1, l2, conv_stop = 32, 1.0, 0.1, 0.1, 255
        tau, gamma, r0, rho = 2.0, 0.01, 0.38, 1.2

        irf, time_axis = model_irf(
            n_channels=n_channels,
            period=period,
            irf_position_p=irf_position_p,
            irf_position_s=irf_position_s,
            irf_width=irf_width
        )
        dt = time_axis[1] - time_axis[0]
        conv_stop = min(len(time_axis) // 2 - 1, conv_stop)
        corrections = np.array([period, g, l1, l2, conv_stop])
        # compute a model function that is later used as "data"
        bg = np.zeros_like(time_axis)
        data = [
            0, 0, 0, 1, 9, 7, 5, 5, 5, 2, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 2, 2, 2, 3, 0, 1, 0,
            1, 1, 1, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        ]
        # create MParam structure that contains all parameters for fitting
        bifl_scatter = 1  # if smaller than zero use soft-bifl fit
        p_2s = 0  # if bigger than zero use p2s_twoIstar
        x = np.zeros(8, dtype=np.float64)
        x[:6] = [tau, gamma, r0, rho, bifl_scatter, p_2s]

        # test targetf23
        m_param = fit2x.CreateMParam(
            irf=irf,  # numpy array, float; instrument response function
            background=bg,  # numpy array, float; background pattern
            data=data,  # numpy array, integer; experimental data / decay
            corrections=corrections,  # numpy array, float; e.g. g-factor
            dt=dt  # float; time interval between micro time bins
        )
        fit2x.compute_signal_and_background(m_param)
        fit2x.targetf23(x, m_param)
        s = np.array([v for v in m_param.get_model()])
        model_ref = np.array(
            [1.55512239e-06, 1.28478411e-06, 1.84127094e-03, 6.97704393e-01
                , 5.26292095e+00, 7.46022080e+00, 5.86203555e+00, 4.37710024e+00
                , 3.27767338e+00, 2.46136154e+00, 1.85332949e+00, 1.39904465e+00
                , 1.05863172e+00, 8.02832214e-01, 6.10103966e-01, 4.64532309e-01
                , 3.54320908e-01, 2.70697733e-01, 2.07119266e-01, 1.58689654e-01
                , 1.21735290e-01, 9.34921381e-02, 7.18751573e-02, 5.53076815e-02
                , 4.25947583e-02, 3.28288183e-02, 2.53192074e-02, 1.95393890e-02
                , 1.50872754e-02, 1.16553438e-02, 9.00806647e-03, 6.96482554e-03
                , 1.32357360e-06, 1.00072013e-05, 2.67318412e-02, 1.32397314e+00
                , 3.09121580e+00, 3.29547010e+00, 2.64779127e+00, 2.11040744e+00
                , 1.67602309e+00, 1.32697773e+00, 1.04788208e+00, 8.25634384e-01
                , 6.49268597e-01, 5.09724532e-01, 3.99592220e-01, 3.12860291e-01
                , 2.44684015e-01, 1.91180006e-01, 1.49249476e-01, 1.16429072e-01
                , 9.07668453e-02, 7.07203018e-02, 5.50733629e-02, 4.28692191e-02
                , 3.33563629e-02, 2.59454278e-02, 2.01748116e-02, 1.56833931e-02
                , 1.21889419e-02, 9.47107726e-03, 7.35784602e-03, 5.71517161e-03]
        )
        # print(s)
        # import pylab as p
        # p.plot(s)
        # p.plot(model_ref)
        # p.show()
        self.assertEqual(
            np.allclose(model_ref, s),
            True
        )

    def test_fit23(self):
        # setup some parameters
        n_channels = 32
        n_photons = 60
        irf_position_p = 2.0
        irf_position_s = 18.0
        irf_width = 0.25
        period, g, l1, l2, conv_stop = 32, 1.0, 0.1, 0.1, 255
        tau, gamma, r0, rho = 2.0, 0.01, 0.38, 1.2
        np.random.seed(0)

        irf_np, time_axis = model_irf(
            n_channels=n_channels,
            period=period,
            irf_position_p=irf_position_p,
            irf_position_s=irf_position_s,
            irf_width=irf_width
        )
        dt = time_axis[1] - time_axis[0]
        conv_stop = min(len(time_axis) // 2 - 1, conv_stop)
        corrections = np.array([period, g, l1, l2, conv_stop])

        # compute a model function that is later used as "data"
        data = [
            0, 0, 0, 1, 9, 7, 5, 5, 5, 2, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 2, 2, 2, 3, 0, 1, 0,
            1, 1, 1, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        ]
        bg = np.zeros_like(irf_np)

        # create MParam structure that contains all parameters for fitting
        m_param = fit2x.CreateMParam(
            irf=irf_np,
            background=bg,
            corrections=corrections,
            dt=dt,
            data=data
        )
        bifl_scatter = -1
        p_2s = 0
        tau = 4.2
        x = np.zeros(8, dtype=np.float64)
        x[:6] = [tau, gamma, r0, rho, bifl_scatter, p_2s]
        fixed = np.array([0, 1, 1, 1], dtype=np.int16)  # lifetime fitted
        twoIstar = fit2x.fit23(x, fixed, m_param)
        fit_res = np.array(
            [1.79364114, 0.01, 0.38, 1.2, -1.,
             0., 0.25974026, 0.25974026]
        )
        # m = np.array([m for m in m_param.get_model()])
        # import pylab as p
        # p.plot(m)
        # p.plot(data)
        # p.plot([x for x in m_param.get_data()])
        # p.show()
        self.assertEqual(np.allclose(fit_res, x), True)
        self.assertAlmostEqual(twoIstar, 0.5309172991292537)

    def test_fit23_2(self):
        irf = np.array(
            [0, 0, 0, 260, 1582, 155, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 22, 1074, 830, 10, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0], dtype=np.float64
        )
        data = np.array(
            [
                0, 0, 0, 1, 9, 7, 5, 5, 5, 2, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 2, 2, 2, 3, 0, 1, 0,
                1, 1, 1, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ]
        )
        settings = {
            'dt': 0.5079365079365079,
            'g_factor': 1.0,
            'l1': 0.1,
            'l2': 0.2,
            'period': 32.0,
            'convolution_stop': 31,
            'irf': irf,
            'background': np.zeros_like(irf)
        }
        fit23 = fit2x.Fit23(**settings)
        tau, gamma, r0, rho = 2.2, 0.01, 0.38, 1.22
        x0 = np.array([tau, gamma, r0, rho])
        fixed = np.array([0, 1, 1, 0])
        r = fit23(
            data=data,
            initial_values=x0,
            fixed=fixed
        )
        data = fit23.data
        model = fit23.model
        self.assertEqual(
            np.allclose(
                r['x'],
                np.array([1.74493538, 0.01, 0.38, 8.75202697, -1.,
                          0., 0.31683168, 0.31683168])
            ),
            True
        )
        self.assertEqual(
            True,
            (data - model).sum() < 0.6
        )

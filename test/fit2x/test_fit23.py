from __future__ import division

import unittest
import numpy as np

import fit2x
from compute_irf import model_irf


class Tests(unittest.TestCase):

    # OK
    def test_modelf23(self):
        # import fit2x
        # import numpy as np
        # import pylab as plt
        period, g, l1, l2, conv_stop = 32., 1.0, 0.1, 0.1, 63
        irf, time_axis = model_irf(
            n_channels=64,
            period=period,
            irf_position_p=2.0,
            irf_position_s=18.0,
            irf_width=0.25
        )
        bg = np.zeros_like(irf) + 0.2
        tau, gamma, r0, rho = 4.0, 0.01, 0.38, 10000.2
        np.random.seed(0)
        param = np.array([tau, gamma, r0, rho])
        conv_stop = min(len(time_axis), conv_stop)
        corrections = np.array([period, g, l1, l2, conv_stop])
        dt = time_axis[1] - time_axis[0]
        out = np.zeros_like(bg)
        fit2x.modelf23(param, irf, bg, dt, corrections, out)
        out_ref = np.array([0.0020231 , 0.00202169, 0.00202036, 0.00201916, 0.00202142,
                            0.00213487, 0.00355841, 0.01042463, 0.02419025, 0.03623856,
                            0.0399039 , 0.03863233, 0.03646924, 0.03436654, 0.03239031,
                            0.03053472, 0.02879244, 0.02715653, 0.02562052, 0.02417828,
                            0.02282411, 0.02155263, 0.02035877, 0.01923782, 0.0181853 ,
                            0.01719705, 0.01626915, 0.01539789, 0.01457984, 0.01381173,
                            0.01309053, 0.01241336, 0.01177754, 0.01118053, 0.01061999,
                            0.01009366, 0.00959948, 0.00913547, 0.00869978, 0.00829071,
                            0.00790661, 0.00754596, 0.00720733, 0.00688938, 0.00659084,
                            0.00631053, 0.00604734, 0.00580021, 0.00556818, 0.00535031,
                            0.00514575, 0.00495367, 0.00477333, 0.00460399, 0.004445  ,
                            0.00429571, 0.00415554, 0.00402392, 0.00390034, 0.00378431,
                            0.00367537, 0.00357307, 0.00347702, 0.00338684, 0.00201   ,
                            0.00200939, 0.00200882, 0.00200846, 0.00201797, 0.00222095,
                            0.00378474, 0.00859881, 0.01498067, 0.01852531, 0.01874061,
                            0.01785824, 0.01689628, 0.0159871 , 0.01513332, 0.01433166,
                            0.01357893, 0.01287214, 0.0122085 , 0.01158537, 0.01100027,
                            0.01045089, 0.00993504, 0.00945068, 0.00899589, 0.00856886,
                            0.00816789, 0.0077914 , 0.00743789, 0.00710596, 0.00679429,
                            0.00650164, 0.00622686, 0.00596885, 0.00572659, 0.00549911,
                            0.00528553, 0.00508498, 0.00489667, 0.00471985, 0.00455383,
                            0.00439794, 0.00425157, 0.00411413, 0.00398509, 0.00386392,
                            0.00375014, 0.00364331, 0.003543  , 0.00344882, 0.00336038,
                            0.00327734, 0.00319937, 0.00312616, 0.00305742, 0.00299288,
                            0.00293227, 0.00287536, 0.00282193, 0.00277176, 0.00272465,
                            0.00268042, 0.00263888, 0.00259989])
        self.assertEqual(
            np.allclose(out, out_ref),
            True
        )

    # OK
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

    # OK
    def test_target23(self):
        n_channels = 32
        irf_position_p = 2.0
        irf_position_s = 18.0
        irf_width = 0.25
        period, g, l1, l2, conv_stop = 32, 1.0, 0.05, 0.05, 128
        tau, gamma, r0, rho = 4.0, 0.01, 0.38, 0.2

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
        model_ref = np.array([1.66519507e-03, 1.46615080e-03, 2.28617854e-03, 3.78693150e-01,
                              2.92975666e+00, 4.66579661e+00, 4.16974650e+00, 3.56147676e+00,
                              3.04670169e+00, 2.61052067e+00, 2.24027400e+00, 1.92543973e+00,
                              1.65725865e+00, 1.42842712e+00, 1.23284490e+00, 1.06540818e+00,
                              9.21839314e-01, 7.98546811e-01, 6.92509823e-01, 6.01182921e-01,
                              5.22417392e-01, 4.54396131e-01, 3.95579697e-01, 3.44661551e-01,
                              3.00530863e-01, 2.62241546e-01, 2.28986446e-01, 2.00075786e-01,
                              1.74919135e-01, 1.53010290e-01, 1.33914598e-01, 1.17258284e-01,
                              1.50297466e-03, 1.32827434e-03, 1.40551369e-02, 6.42147385e-01,
                              1.62959724e+00, 1.92890467e+00, 1.75911255e+00, 1.59382873e+00,
                              1.43975395e+00, 1.29717655e+00, 1.16603916e+00, 1.04604027e+00,
                              9.36711296e-01, 8.37475283e-01, 7.47691217e-01, 6.66687042e-01,
                              5.93783899e-01, 5.28313546e-01, 4.69630533e-01, 4.17120352e-01,
                              3.70204542e-01, 3.28343503e-01, 2.91037617e-01, 2.57827136e-01,
                              2.28291191e-01, 2.02046204e-01, 1.78743905e-01, 1.58069115e-01,
                              1.39737410e-01, 1.23492763e-01, 1.09105210e-01, 9.63685926e-02])
        # import pylab as p
        # p.semilogy(s, label='Model')
        # p.semilogy(model_ref, label='Reference')
        # p.legend()
        # p.show()
        self.assertEqual(
            np.allclose(model_ref, s),
            True
        )

    def test_fit23(self):
        # setup some parameters
        n_channels = 32
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
        fit_res = np.array([ 1.79364115,  0.01      ,  0.38      ,  1.2       , -1.        ,
                             0.        ,  0.25974026,  0.25974026])
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

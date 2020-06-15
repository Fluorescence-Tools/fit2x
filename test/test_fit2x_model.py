from __future__ import division

import unittest
import numpy as np
import scipy.stats

import fit2x


def model_irf(
        n_channels: int = 256,
        period: float = 32,
        irf_position_p: float = 2.0,
        irf_position_s: float = 18.0,
        irf_width: float = 0.25
):
    time_axis = np.linspace(0, period, n_channels * 2)
    irf_np = scipy.stats.norm.pdf(time_axis, loc=irf_position_p,
                                  scale=irf_width) + \
             scipy.stats.norm.pdf(time_axis, loc=irf_position_s,
                                  scale=irf_width)
    return irf_np, time_axis


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

    # def test_modelf24(self):
    #     irf, time_axis = model_irf(
    #         n_channels=64,
    #         period=32.,
    #         irf_position_p=2.0,
    #         irf_position_s=18.0,
    #         irf_width=0.25
    #     )
    #     bg = np.zeros_like(irf) + 0.2
    #
    #     tau1 = 4.0
    #     gamma = 0.01
    #     tau2 = 0.5
    #     a2 = 0.9
    #     offset = 1
    #     param = np.array([tau1, gamma, tau2, a2, offset])
    #     period, g, l1, l2, conv_stop = 32, 1.0, 0.1, 0.1, 63
    #     conv_stop = min(len(time_axis), conv_stop)
    #     corrections = np.array([period, g, l1, l2, conv_stop])
    #     dt = time_axis[1] - time_axis[0]
    #     out = np.zeros_like(bg)
    #     fit2x.modelf24(param, irf, bg, dt, corrections, out)
    #
    #     out_ref = np.array(
    #         [0.01571078, 0.01571031, 0.01570988, 0.01570959, 0.01572048,
    #          0.01609258, 0.02058282, 0.04054985, 0.07322865, 0.08840347,
    #          0.07745669, 0.05963912, 0.04647572, 0.03811999, 0.0328345,
    #          0.02942218, 0.02715495, 0.02559208, 0.02446654, 0.02361628,
    #          0.02294272, 0.02238568, 0.0219082, 0.02148738, 0.02110889,
    #          0.02076354, 0.02044534, 0.02015024, 0.01987537, 0.01961861,
    #          0.01937834, 0.01915323, 0.01894215, 0.01874414, 0.01855832,
    #          0.01838392, 0.0182202, 0.0180665, 0.01792219, 0.01778671,
    #          0.0176595, 0.01754006, 0.01742791, 0.01732261, 0.01722374,
    #          0.01713091, 0.01704375, 0.0169619, 0.01688505, 0.0168129,
    #          0.01674515, 0.01668153, 0.0166218, 0.01656572, 0.01651306,
    #          0.01646361, 0.01641718, 0.01637359, 0.01633266, 0.01629423,
    #          0.01625814, 0.01622426, 0.01619244, 0.01616257, 0.01571054,
    #          0.01571009, 0.01570967, 0.01571062, 0.01578379, 0.01724953,
    #          0.02796797, 0.05699103, 0.08457921, 0.08515608, 0.06828081,
    #          0.05234779, 0.04181477, 0.03518192, 0.03094818, 0.0281781,
    #          0.02630526, 0.02498664, 0.02401429, 0.02326186, 0.02265236,
    #          0.02213868, 0.02169176, 0.02129351, 0.0209325, 0.02060134,
    #          0.02029511, 0.02001042, 0.01974484, 0.01949651, 0.01926396,
    #          0.019046, 0.01884157, 0.01864976, 0.01846974, 0.01830077,
    #          0.01814213, 0.01799321, 0.01785338, 0.0177221, 0.01759884,
    #          0.0174831, 0.01737443, 0.0172724, 0.0171766, 0.01708664,
    #          0.01700218, 0.01692287, 0.01684841, 0.01677849, 0.01671284,
    #          0.0166512, 0.01659332, 0.01653897, 0.01648795, 0.01644003,
    #          0.01639505, 0.0163528, 0.01631314, 0.0162759, 0.01624093,
    #          0.0162081, 0.01617727, 0.01614833]
    #     )
    #     self.assertEqual(
    #         np.allclose(out, out_ref),
    #         True
    #     )

    def test_LVDoubleArray(self):
        lv_array = fit2x.CreateLVDoubleArray(10)
        self.assertEqual(
            len(lv_array), 10
        )
        lv_array[0] = 10
        vec = np.linspace(1, 10., 10, dtype=np.float)
        for i, v in enumerate(vec):
            lv_array[i] = v
        l = [v for v in lv_array]
        self.assertListEqual(
            list(lv_array), l
        )

    @unittest.expectedFailure
    def test_LVDoubleArray_slice(self):
        lv_array = fit2x.CreateLVDoubleArray(10)
        vec = np.linspace(1, 10., 10, dtype=np.float)
        # will fail slicing not supported
        lv_array[:] = vec

    def test_LVI32Array(self):
        lv_array = fit2x.CreateLVI32Array(10)
        self.assertEqual(
            len(lv_array), 10
        )
        lv_array[0] = 10
        vec = np.linspace(1, 10, 10, dtype=np.int)
        # will fail slicing not supported
        # lv_array[:] = vec
        for i, v in enumerate(vec):
            lv_array[i] = int(v)
        l = [v for v in lv_array]
        self.assertListEqual(
            list(lv_array), l
        )

    @unittest.expectedFailure
    def test_LVI32Array_slice(self):
        lv_array = fit2x.CreateLVI32Array(10)
        vec = np.linspace(1, 10, 10, dtype=np.int)
        lv_array[:] = vec

    # def test_MParam_1(self):
    #     # create empty MParam structure
    #     n_channels = 256
    #     n_corrections = 4
    #     dt = 0.032
    #     parameter_group = fit2x.CreateMParam(
    #         n_channels=n_channels,
    #         n_corrections=n_corrections,
    #         dt=dt
    #     )
    #     irf = parameter_group.get_irf()
    #     bg = parameter_group.get_background()
    #     corrections = parameter_group.get_corrections()
    #     self.assertEqual(len(irf), 2 * n_channels)
    #     self.assertEqual(len(bg), 2 * n_channels)
    #     self.assertEqual(type(irf), fit2x.LVDoubleArray)
    #     self.assertEqual(type(bg), fit2x.LVDoubleArray)
    #     self.assertEqual(type(parameter_group.get_model()), fit2x.LVDoubleArray)
    #     self.assertEqual(parameter_group.dt, dt)
    #     irf[1] = 22.0
    #     irf2 = parameter_group.get_irf()
    #     self.assertEqual(irf2[1], 22.0)

    def test_MParam_2(self):
        # create filled MParam structure
        n_corrections = 5
        corrections_np = np.zeros(n_corrections, dtype=np.float)
        irf_np, timeaxis = model_irf()
        dt = timeaxis[1] - timeaxis[0]
        bg_np = np.zeros_like(irf_np)
        parameter_group = fit2x.CreateMParam(
            irf=irf_np,
            background=bg_np,
            corrections=corrections_np,
            dt=dt,
        )
        irf = parameter_group.get_irf()
        bg = parameter_group.get_background()
        corrections = parameter_group.get_corrections()
        self.assertEqual(len(irf_np), len(irf))
        self.assertEqual(len(bg_np), len(bg))
        irf[0] = 123
        self.assertEqual(irf[0], 123)
        bg[0] = 123
        self.assertEqual(bg[0], 123)
        corrections[0] = 11
        self.assertEqual(corrections[0], 11)

    def test_lv_param(self):
        for i in range(200):
            lv_array = fit2x.CreateLVDoubleArray(10)
            a = np.ones(111, dtype=np.float)
            lv_array.set_data(a)
            x = [x for x in lv_array]
            self.assertListEqual(list(a), x)

        for i in range(200):
            lv_array = fit2x.CreateLVI32Array(10)
            a = np.ones(111, dtype=np.int32)
            lv_array.set_data(a)
            x = [x for x in lv_array]
            self.assertListEqual(list(a), x)

        for i in range(1000):
            dt = 0.032
            a = np.ones(111, dtype=np.float)
            parameter_group = fit2x.CreateMParam(
                irf=a,
                background=a,
                corrections=a,
                dt=dt
            )
            irf = parameter_group.get_irf()
            irf.set_data(a)
            x = [x for x in irf]
            self.assertListEqual(list(a), x)

            bg = parameter_group.get_background()
            bg.set_data(a)
            x = [x for x in bg]
            self.assertListEqual(list(a), x)

            corrections = parameter_group.get_corrections()
            corrections.set_data(a)
            x = [x for x in corrections]
            self.assertListEqual(list(a), x)

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
        conv_stop = len(time_axis) / 2
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
            irf=irf_np,
            background=bg,
            data=data,
            corrections=corrections,
            dt=dt
        )
        fit2x.compute_signal_and_background(m_param)
        fit2x.targetf23(x, m_param)
        s = np.array([v for v in m_param.get_model()])
        model_ref = np.array(
            [1.63983931e-06, 1.35489900e-06, 1.94264755e-03, 7.35846365e-01,
             5.44784341e+00, 7.10922785e+00, 5.24679453e+00, 3.76146623e+00,
             2.76080690e+00, 2.06183293e+00, 1.55887332e+00, 1.18865676e+00,
             9.11602173e-01, 7.01831779e-01, 5.41723368e-01, 4.18852552e-01,
             3.24214237e-01, 2.51144444e-01, 1.94637102e-01, 1.50891873e-01,
             1.17002911e-01, 9.07375112e-02, 7.03746048e-02, 5.45846545e-02,
             4.23391231e-02, 3.28415836e-02, 2.54749584e-02, 1.97609386e-02,
             1.53286789e-02, 1.18906035e-02, 9.22368335e-03, 7.15493571e-03,
             1.39676185e-06, 1.05584561e-05, 2.82050670e-02, 1.40058981e+00,
             3.44379723e+00, 3.83150593e+00, 3.15607116e+00, 2.53652297e+00,
             2.01248245e+00, 1.58390511e+00, 1.24023852e+00, 9.67954820e-01,
             7.53846292e-01, 5.86288813e-01, 4.55564955e-01, 3.53780944e-01,
             2.74632691e-01, 2.13138185e-01, 1.65386149e-01, 1.28318843e-01,
             9.95522889e-02, 7.72310792e-02, 5.99128375e-02, 4.64771100e-02,
             3.60539405e-02, 2.79680824e-02, 2.16955288e-02, 1.68296938e-02,
             1.30551288e-02, 1.01271077e-02, 7.85577899e-03, 6.09386452e-03]
        )
        # import pylab as p
        # p.plot(s)
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
        param = np.array([tau, gamma, r0, rho])
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
        x = np.zeros(8, dtype=np.float64)
        x[:6] = [4.2, gamma, r0, rho, bifl_scatter, p_2s]

        # test fitting, lifetime free
        fixed = np.array([0, 1, 1, 1], dtype=np.int16)
        Istar = fit2x.fit23(x, fixed, m_param)
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
        self.assertAlmostEqual(Istar, 0.5309172991292537)


from __future__ import division

import unittest
import numpy as np
import scipy

import fit2x
from compute_irf import model_irf


class Tests(unittest.TestCase):

    data = np.array(
        [
            0, 0, 0, 1, 9, 7, 5, 5, 5, 2, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 2, 2, 2, 3, 0, 1, 0,
            1, 1, 1, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        ],
        dtype=np.float
    )
    w_sq = np.sqrt(data)**2.
    model = np.array(
        [
            1.55512239e-06, 1.28478411e-06, 1.84127094e-03, 6.97704393e-01,
            5.26292095e+00, 7.46022080e+00, 5.86203555e+00, 4.37710024e+00,
            3.27767338e+00, 2.46136154e+00, 1.85332949e+00, 1.39904465e+00,
            1.05863172e+00, 8.02832214e-01, 6.10103966e-01, 4.64532309e-01,
            3.54320908e-01, 2.70697733e-01, 2.07119266e-01, 1.58689654e-01,
            1.21735290e-01, 9.34921381e-02, 7.18751573e-02, 5.53076815e-02,
            4.25947583e-02, 3.28288183e-02, 2.53192074e-02, 1.95393890e-02,
            1.50872754e-02, 1.16553438e-02, 9.00806647e-03, 6.96482554e-03,
            1.32357360e-06, 1.00072013e-05, 2.67318412e-02, 1.32397314e+00,
            3.09121580e+00, 3.29547010e+00, 2.64779127e+00, 2.11040744e+00,
            1.67602309e+00, 1.32697773e+00, 1.04788208e+00, 8.25634384e-01,
            6.49268597e-01, 5.09724532e-01, 3.99592220e-01, 3.12860291e-01,
            2.44684015e-01, 1.91180006e-01, 1.49249476e-01, 1.16429072e-01,
            9.07668453e-02, 7.07203018e-02, 5.50733629e-02, 4.28692191e-02,
            3.33563629e-02, 2.59454278e-02, 2.01748116e-02, 1.56833931e-02,
            1.21889419e-02, 9.47107726e-03, 7.35784602e-03, 5.71517161e-03
        ]
    )

    def test_rescale(self):
        model = np.copy(self.model)
        data = np.copy(self.data)
        model /= 10
        scale = fit2x.rescale(
            fit=model,
            decay=data,
            start=0,
            stop=-1
        )
        print(scale)
        self.assertAlmostEqual(
            scale, 10.1010, places=2
        )

    def test_rescale_w(self):
        model = np.copy(self.model)
        data = np.copy(self.data)
        w_sq = np.copy(self.w_sq)
        model /= 10
        scale = fit2x.rescale_w(
            fit=model,
            decay=data,
            w_sq=w_sq,
            start=0,
            stop=-1
        )
        self.assertAlmostEqual(
            scale, 10, places=1
        )

    def test_rescale_w_bg(self):
        model = np.copy(self.model)
        data = np.copy(self.data)
        w_sq = np.copy(self.w_sq)
        model /= 10
        scale = fit2x.rescale_w_bg(
            fit=model,
            decay=data,
            w_sq=w_sq,
            bg=1.0,
            start=0,
            stop=-1
        )
        self.assertAlmostEqual(
            scale, 6.524429195170, places=2
        )

    def test_fconv(self):
        period = 12.0
        lifetime_spectrum = np.array([1.0, 4.1])
        irf, time_axis = model_irf(
            n_channels=32,
            period=period,
            irf_position_p=2.0,
            irf_position_s=2.0,
            irf_width=0.25
        )
        dt = time_axis[1] - time_axis[0]

        irf += 0.05
        model_ref = np.array([0.        , 0.00930764, 0.01819277, 0.02667457, 0.03477278,
                              0.04255003, 0.05082231, 0.06731845, 0.13140084, 0.34023069,
                              0.76778777, 1.2947723 , 1.67334688, 1.80724789, 1.79049142,
                              1.72765984, 1.65940007, 1.59342535, 1.53040088, 1.47023605,
                              1.41280237, 1.35797588, 1.30563825, 1.25567648, 1.20798272,
                              1.16245403, 1.1189921 , 1.07750313, 1.03789755, 1.00008986,
                              0.96399845, 0.92954541, 0.89665637, 0.86526032, 0.83528949,
                              0.80667919, 0.77936766, 0.75329593, 0.72840772, 0.70464932,
                              0.68196942, 0.66031908, 0.63965156, 0.61992224, 0.60108853,
                              0.58310977, 0.56594717, 0.54956365, 0.53392387, 0.51899405,
                              0.50474198, 0.49113687, 0.47814937, 0.46575143, 0.4539163 ,
                              0.44261843, 0.43183342, 0.421538  , 0.41170994, 0.40232802,
                              0.39337199, 0.38482252, 0.37666116, 0.36887028])
        model_fconv = np.zeros_like(irf)
        fit2x.fconv(
            fit=model_fconv,
            irf=irf,
            x=lifetime_spectrum,
            dt=dt
        )
        self.assertEqual(
            np.allclose(model_ref, model_fconv), True
        )

    def test_fconv_per(self):
        period = 12.0
        lifetime_spectrum = np.array([1.0, 4.1])
        irf, time_axis = model_irf(
            n_channels=32,
            period=period,
            irf_position_p=2.0,
            irf_position_s=2.0,
            irf_width=0.25
        )
        dt = time_axis[1] - time_axis[0]

        model_fconv_per = np.zeros_like(irf)
        fit2x.fconv_per(
            fit=model_fconv_per,
            irf=irf,
            x=lifetime_spectrum,
            period=12,
            start=5,
            stop=-1,
            dt=dt
        )
        ref = np.array([0.00000000e+00, 1.28147164e-12, 2.39814805e-10, 2.51323110e-08,
                        1.48317225e-06, 1.39826223e-01, 1.34375008e-01, 1.37770632e-01,
                        1.89347207e-01, 3.86238943e-01, 8.02399833e-01, 1.31850551e+00,
                        1.68669508e+00, 1.81068251e+00, 1.78446248e+00, 1.71259694e+00,
                        1.63571330e+00, 1.56150620e+00, 1.49062306e+00, 1.42295629e+00,
                        1.35836123e+00, 1.29669846e+00, 1.23783487e+00, 1.18164339e+00,
                        1.12800272e+00, 1.07679707e+00, 1.02791590e+00, 9.81253684e-01,
                        9.36709701e-01, 8.94187790e-01, 8.53596160e-01, 8.14847186e-01,
                        7.77857221e-01, 7.42546414e-01, 7.08838540e-01, 6.76660834e-01,
                        6.45943835e-01, 6.16621232e-01, 5.88629729e-01, 5.61908899e-01,
                        5.36401060e-01, 5.12051149e-01, 4.88806601e-01, 4.66617239e-01,
                        4.45435162e-01, 4.25214646e-01, 4.05912039e-01, 3.87485673e-01,
                        3.69895772e-01, 3.53104363e-01, 3.37075201e-01, 3.21773681e-01,
                        3.07166774e-01, 2.93222947e-01, 2.79912099e-01, 2.67205497e-01,
                        2.55075711e-01, 2.43496557e-01, 2.32443037e-01, 2.21891292e-01,
                        2.11818543e-01, 2.02203046e-01, 1.93024044e-01, 1.84261723e-01])
        self.assertEqual(
            np.allclose(model_fconv_per, ref), True
        )

    def test_fconv_per_cs(self):
        period = 12.0
        lifetime_spectrum = np.array([1.0, 4.1])
        irf, time_axis = model_irf(
            n_channels=32,
            period=period,
            irf_position_p=2.0,
            irf_position_s=2.0,
            irf_width=0.25
        )
        dt = time_axis[1] - time_axis[0]
        model_fconv_per_cs = np.zeros_like(irf)
        fit2x.fconv_per_cs(
            fit=model_fconv_per_cs,
            irf=irf,
            x=lifetime_spectrum,
            period=12,
            stop=25,
            conv_stop=50,
            dt=dt
        )
        ref = np.array([0.03336827, 0.03296381, 0.03256425, 0.03216955, 0.03178   ,
                        0.03140744, 0.0312624 , 0.0333835 , 0.04823524, 0.10202368,
                        0.2160135 , 0.36090352, 0.47135935, 0.51967336, 0.52850075,
                        0.52456987, 0.51844462, 0.51217302, 0.50596532, 0.49983249,
                        0.49377399, 0.48778892, 0.4818764 , 0.47603555, 0.47026549,
                        0.46456537, 0.43463299, 0.42936477, 0.42416042, 0.41901914,
                        0.41394018, 0.40892279, 0.40396621, 0.39906971, 0.39423256,
                        0.38945404, 0.38473345, 0.38007007, 0.37546322, 0.3709122 ,
                        0.36641635, 0.361975  , 0.35758748, 0.35325314, 0.34897133,
                        0.34474143, 0.3405628 , 0.33643481, 0.33235687, 0.32832835,
                        0.32434866, 0.32041721, 0.31653341, 0.31269669, 0.30890647,
                        0.3051622 , 0.30146331, 0.29780925, 0.29419949, 0.29063348,
                        0.28711069, 0.2836306 , 0.2801927 , 0.27679647])
        self.assertEqual(
            np.allclose(ref, model_fconv_per_cs), True
        )

    def test_sconv(self):
        period = 12.0
        irf, time_axis = model_irf(
            n_channels=32,
            period=period,
            irf_position_p=2.0,
            irf_position_s=2.0,
            irf_width=0.25
        )
        tau = 4.1
        decay = np.exp(-time_axis / tau)
        model_sconv = np.zeros_like(irf)
        fit2x.sconv(
            fit=model_sconv,
            irf=irf,
            model=decay
        )
        ref = np.array([0.00000000e+00, 6.72772613e-12, 1.25902772e-09, 1.31944633e-07,
                        7.78665430e-06, 2.60011255e-04, 4.95321954e-03, 5.45801295e-02,
                        3.55713481e-01, 1.41837340e+00, 3.63088091e+00, 6.36684281e+00,
                        8.32504638e+00, 9.00004436e+00, 8.88536083e+00, 8.52999555e+00,
                        8.14728984e+00, 7.77768563e+00, 7.42462513e+00, 7.08758461e+00,
                        6.76584390e+00, 6.45870861e+00, 6.16551571e+00, 5.88563229e+00,
                        5.61845418e+00, 5.36340461e+00, 5.11993301e+00, 4.88751379e+00,
                        4.66564524e+00, 4.45384840e+00, 4.25166608e+00, 4.05866182e+00,
                        3.87441898e+00, 3.69853983e+00, 3.53064472e+00, 3.37037120e+00,
                        3.21737330e+00, 3.07132073e+00, 2.93189821e+00, 2.79880477e+00,
                        2.67175310e+00, 2.55046895e+00, 2.43469048e+00, 2.32416777e+00,
                        2.21866224e+00, 2.11794613e+00, 2.02180203e+00, 1.93002238e+00,
                        1.84240907e+00, 1.75877296e+00, 1.67893351e+00, 1.60271837e+00,
                        1.52996301e+00, 1.46051039e+00, 1.39421056e+00, 1.33092041e+00,
                        1.27050331e+00, 1.21282885e+00, 1.15777252e+00, 1.10521547e+00,
                        1.05504424e+00, 1.00715054e+00, 9.61430968e-01, 9.17786836e-01])
        self.assertEqual(
            np.allclose(ref, model_sconv), True
        )

    def test_convolve_lifetime_spectrum_periodic(self):
        time_axis = np.linspace(0, 25, 25)
        irf_position = 6.0
        irf_width = 1.0
        irf = scipy.stats.norm.pdf(time_axis, loc=irf_position, scale=irf_width)
        lifetime_spectrum = np.array([0.2, 1.1, 0.8, 4.0])
        model_decay = np.zeros_like(time_axis)
        fit2x.fconv_per_cs_time_axis(
            model=model_decay,
            time_axis=time_axis,
            lifetime_spectrum=lifetime_spectrum,
            instrument_response_function=irf,
            period=16.0
        )
        reference = np.array(
            [0.00560653, 0.00432208, 0.00342868, 0.00603476, 0.0454072,
             0.21058509, 0.4551221, 0.55543513, 0.48165047, 0.36885986,
             0.27946424, 0.21331303, 0.16359693, 0.12577494, 0.09681668,
             0.07457228, 0.05745678, 0.04427657, 0.03412254, 0.02629821,
             0.02026841, 0.01562132, 0.01203976, 0.00927939, 0.0071519])
        self.assertEqual(
            np.allclose(reference, model_decay), True
        )


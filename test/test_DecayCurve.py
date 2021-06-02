from __future__ import division
import unittest

import numpy as np
import scipy.stats
import fit2x
try:
    import tttrlib
    skip_tttr_tests = False
except ImportError:
    skip_tttr_tests = True


class Tests(unittest.TestCase):

    def test_shift(self):
        time_axis = np.linspace(0, 12, 25)
        irf_position = 6.0
        irf_width = 1.0
        irf = scipy.stats.norm.pdf(time_axis, loc=irf_position, scale=irf_width)
        # integer shift
        shifted_irf_ip = fit2x.DecayCurve.shift_array(
            input=irf,
            shift=1.0
        )
        ref = np.array(
            [0.00000000e+00, 1.48671951e-06, 1.59837411e-05, 1.33830226e-04,
             8.72682695e-04, 4.43184841e-03, 1.75283005e-02, 5.39909665e-02,
             1.29517596e-01, 2.41970725e-01, 3.52065327e-01, 3.98942280e-01,
             3.52065327e-01, 2.41970725e-01, 1.29517596e-01, 5.39909665e-02,
             1.75283005e-02, 4.43184841e-03, 8.72682695e-04, 1.33830226e-04,
             1.59837411e-05, 1.48671951e-06, 1.07697600e-07, 6.07588285e-09,
             6.07588285e-09]
        )
        self.assertEqual(np.allclose(shifted_irf_ip, ref), True)
        shifted_irf_in = fit2x.DecayCurve.shift_array(
            input=irf,
            shift=-1.0
        )
        ref = np.array(
            [6.07588285e-09, 6.07588285e-09, 1.07697600e-07, 1.48671951e-06,
             1.59837411e-05, 1.33830226e-04, 8.72682695e-04, 4.43184841e-03,
             1.75283005e-02, 5.39909665e-02, 1.29517596e-01, 2.41970725e-01,
             3.52065327e-01, 3.98942280e-01, 3.52065327e-01, 2.41970725e-01,
             1.29517596e-01, 5.39909665e-02, 1.75283005e-02, 4.43184841e-03,
             8.72682695e-04, 1.33830226e-04, 1.59837411e-05, 1.48671951e-06,
             1.07697600e-07]
        )
        self.assertEqual(np.allclose(shifted_irf_in, ref), True)

        # floating shift
        shifted_irf_fp = fit2x.DecayCurve.shift_array(
            input=irf,
            shift=1.5
        )
        ref = np.array(
            [0.00000000e+00, 8.73523031e-06, 7.49069834e-05, 5.03256460e-04,
             2.65226555e-03, 1.09800745e-02, 3.57596335e-02, 9.17542811e-02,
             1.85744160e-01, 2.97018026e-01, 3.75503804e-01, 3.75503804e-01,
             2.97018026e-01, 1.85744160e-01, 9.17542811e-02, 3.57596335e-02,
             1.09800745e-02, 2.65226555e-03, 5.03256460e-04, 7.49069834e-05,
             8.73523031e-06, 7.97208558e-07, 5.68867416e-08, 6.07588285e-09,
             5.68867416e-08]
        )
        self.assertEqual(np.allclose(ref, shifted_irf_fp), True)

        shifted_irf_fn = fit2x.DecayCurve.shift_array(
            input=irf,
            shift=-1.5
        )
        ref = np.array(
            [6.07588285e-09, 5.68867416e-08, 7.97208558e-07, 8.73523031e-06,
             7.49069834e-05, 5.03256460e-04, 2.65226555e-03, 1.09800745e-02,
             3.57596335e-02, 9.17542811e-02, 1.85744160e-01, 2.97018026e-01,
             3.75503804e-01, 3.75503804e-01, 2.97018026e-01, 1.85744160e-01,
             9.17542811e-02, 3.57596335e-02, 1.09800745e-02, 2.65226555e-03,
             5.03256460e-04, 7.49069834e-05, 8.73523031e-06, 7.97208558e-07,
             5.68867416e-08]
        )
        self.assertEqual(np.allclose(ref, shifted_irf_fn), True)

        # rollover
        shifted_irf_irp = fit2x.DecayCurve.shift_array(
            input=irf,
            shift=26.0
        )
        self.assertEqual(
            np.allclose(shifted_irf_irp, shifted_irf_ip), True
        )

    def test_shift_2(self):
        x = np.arange(0, 8, dtype=np.float64)
        y = np.zeros_like(x)
        y[5] = 1.0
        dc = fit2x.DecayCurve(x, y)
        self.assertEqual(dc.shift, 0.0)
        dc.shift = 2
        self.assertEqual(np.allclose(np.roll(y, -2), dc.y), True)
        dc.shift = 2.5
        self.assertEqual(
            np.allclose(
                np.array([0., 0., 0.5, 0.5, 0., 0., 0., 0.]),
                dc.y), True
        )

    def test_add_irf(self):
        time_axis = np.linspace(0, 10, 64)
        irf_position = 1.0
        irf_width = 0.5
        irf = scipy.stats.norm.pdf(
            time_axis,
            loc=irf_position,
            scale=irf_width
        )
        lifetime_spectrum = np.array([1.0, 4.0])
        model_decay = np.zeros_like(time_axis)
        fit2x.fconv_per_cs_time_axis(
            model_decay,
            time_axis=time_axis,
            lifetime_spectrum=lifetime_spectrum,
            instrument_response_function=irf
        )
        model_incl_irf = fit2x.DecayCurve.add_arrays(
            curve1=model_decay,
            curve2=irf,
            start=0,
            stop=-1,
            areal_fraction_curve2=0.9
        )
        ref = np.array([0.3115402 , 0.55977279, 0.91044342, 1.33818461, 1.778436  ,
                        2.13790286, 2.3256597 , 2.2907235 , 2.04502752, 1.65767663,
                        1.22415371, 0.82903448, 0.5217272 , 0.31314788, 0.18786564,
                        0.12042988, 0.08730147, 0.07188893, 0.06454703, 0.06049367,
                        0.05767012, 0.05529761, 0.05311444, 0.05104094, 0.04905375,
                        0.04714504, 0.04531079, 0.04354795, 0.04185369, 0.04022535,
                        0.03866037, 0.03715627, 0.03571068, 0.03432134, 0.03298605,
                        0.03170271, 0.0304693 , 0.02928388, 0.02814458, 0.0270496 ,
                        0.02599722, 0.02498578, 0.0240137 , 0.02307943, 0.02218152,
                        0.02131853, 0.02048913, 0.01969199, 0.01892586, 0.01818954,
                        0.01748187, 0.01680173, 0.01614805, 0.0155198 , 0.01491599,
                        0.01433568, 0.01377794, 0.0132419 , 0.01272672, 0.01223158,
                        0.01175571, 0.01129834, 0.01085878, 0.01043631])
        self.assertEqual(np.allclose(ref, model_incl_irf), True)

    @unittest.skipIf(skip_tttr_tests, "Cloud import tttrlib")
    def test_DecayCurve_set_tttr(self):
        tttr_fn = "../tttr-data/bh/bh_spc132_sm_dna/m000.spc"
        tttr = tttrlib.TTTR(tttr_fn, "SPC-130")
        curve = fit2x.DecayCurve()
        curve.set_tttr(tttr, 256)
        # plt.semilogy(curve.x * 1e9, curve.y)
        # plt.show()
        self.assertEqual(
            np.allclose(
                curve.y,
                np.array([0., 547., 34014., 74748., 15763., 9930., 7289., 5942.,
                          4909., 4386., 4029., 3683., 3419., 4021., 1758., 0.])

            ),
            True
        )
        self.assertEqual(
            np.allclose(
                curve.x,
                np.array([0.00000000e+00, 1.28746033e-14, 2.57492065e-14, 3.86238098e-14,
                          5.14984131e-14, 6.43730164e-14, 7.72476196e-14, 9.01222229e-14,
                          1.02996826e-13, 1.15871429e-13, 1.28746033e-13, 1.41620636e-13,
                          1.54495239e-13, 1.67369843e-13, 1.80244446e-13, 1.93119049e-13])
            ),
            True
        )

    def test_DecayCurve_init(self):
        x = np.linspace(0, 10, 10)
        y = np.sin(x) * 1000
        curve1 = fit2x.DecayCurve(x, y)
        curve2 = fit2x.DecayCurve(x=x, y=y)
        self.assertEqual(np.allclose(curve1.x, curve2.x), True)
        self.assertEqual(np.allclose(curve1.y, curve2.y), True)
        # Noise is by default Poisson noise
        ey = np.clip(np.sqrt(np.abs(y)), 1, 1e12)
        curve1 = fit2x.DecayCurve(x, y)
        self.assertEqual(np.allclose(curve1.ey, ey), True)

    def test_DecayCurve_operator(self):
        x = np.linspace(0, 10, 10)
        y = np.sin(x)
        ey = np.clip(np.sqrt(np.abs(y)), 1, 1e12)

        # Multiplication
        curve1 = fit2x.DecayCurve(x, y)
        curve2 = curve1 * curve1
        curve3 = curve1 * 2
        curve1 *= 2
        self.assertEqual(np.allclose(curve1.y, y * 2), True)
        self.assertEqual(np.allclose(curve1.ey, ey * 2), True)
        self.assertEqual(np.allclose(curve3.y, curve1.y), True)
        self.assertEqual(np.allclose(curve2.x, curve1.x), True)
        self.assertEqual(np.allclose(curve3.ey, curve1.ey), True)
        self.assertEqual(np.allclose(curve2.y, y * y), True)

        # Addition
        curve1 = fit2x.DecayCurve(x, y)
        curve3 = curve1 + curve1
        curve2 = curve1 + 2
        curve1 += 2
        self.assertEqual(np.allclose(curve1.y, y + 2), True)
        self.assertEqual(np.allclose(curve1.ey, ey), True)
        self.assertEqual(np.allclose(curve2.y, curve1.y), True)
        self.assertEqual(np.allclose(curve2.x, curve1.x), True)
        self.assertEqual(np.allclose(curve2.ey, curve1.ey), True)
        self.assertEqual(np.allclose(curve3.y, y + y), True)

    def test_decay_curve_resize(self):
        d = fit2x.DecayCurve()
        d.x = np.array([0, 0.1])
        d.resize(5)
        self.assertEqual(np.allclose(d.x, np.array([0. , 0.1, 0.2, 0.3, 0.4])), True)
        d.resize(6)
        self.assertEqual(np.allclose(d.x, np.array([0. , 0.1, 0.2, 0.3, 0.4, 0.5])), True)

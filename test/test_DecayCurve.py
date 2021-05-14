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

    # Seems ok
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
        ref = np.array([0.34899169, 0.62778997, 1.02023263, 1.49880469, 1.99114485,
                        2.39274598, 2.60181415, 2.56135465, 2.28487964, 1.84993179,
                        1.36355308, 0.92050288, 0.576091  , 0.34247794, 0.20230614,
                        0.12700672, 0.09017091, 0.07318947, 0.06524602, 0.06097861,
                        0.05807925, 0.05567522, 0.05347352, 0.05138519, 0.04938444,
                        0.04746283, 0.04561622, 0.04384149, 0.04213581, 0.0404965 ,
                        0.03892096, 0.03740672, 0.03595139, 0.03455269, 0.0332084 ,
                        0.03191641, 0.03067468, 0.02948127, 0.02833429, 0.02723193,
                        0.02617246, 0.0251542 , 0.02417557, 0.023235  , 0.02233103,
                        0.02146223, 0.02062723, 0.01982472, 0.01905343, 0.01831215,
                        0.0175997 , 0.01691498, 0.01625689, 0.01562441, 0.01501654,
                        0.01443231, 0.01387081, 0.01333116, 0.01281251, 0.01231403,
                        0.01183495, 0.0113745 , 0.01093197, 0.01050666])
        self.assertEqual(np.allclose(ref, model_incl_irf), True)

    @unittest.skipIf(skip_tttr_tests, "Cloud import tttrlib")
    def test_DecayCurve_set_tttr(self):
        tttr_fn = "../tttr-data/bh/bh_spc132_sm_dna/m000.spc"
        tttr = tttrlib.TTTR(tttr_fn, "SPC-130")

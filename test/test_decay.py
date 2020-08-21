from __future__ import division
import unittest

import scipy.stats
# import tttrlib
import fit2x
import numpy as np

print("Test: ", __file__)


class Tests(unittest.TestCase):

    def test_getter_setter(self):
        decay = fit2x.Decay()

        decay.use_amplitude_threshold = True
        self.assertEqual(decay.use_amplitude_threshold, True)
        decay.use_amplitude_threshold = False
        self.assertEqual(decay.use_amplitude_threshold, False)

        decay.amplitude_threshold = 11
        self.assertEqual(decay.amplitude_threshold, 11)
        decay.amplitude_threshold = 2.2
        self.assertEqual(decay.amplitude_threshold, 2.2)

        decay.constant_offset = 11
        self.assertEqual(decay.constant_offset, 11)
        decay.constant_offset = 2.2
        self.assertEqual(decay.constant_offset, 2.2)

        # This does not work because there is no IRF
        # decay.set_irf_shift_channels(11)
        # self.assertEqual(decay.get_irf_shift_channels(), 11)
        # decay.set_irf_shift_channels(2.2)
        # self.assertEqual(decay.get_irf_shift_channels(), 2.2)

        decay.number_of_photons = 11
        self.assertEqual(decay.number_of_photons, 11)
        decay.number_of_photons = 2.2
        self.assertEqual(decay.number_of_photons, 2.2)

        decay.scatter_fraction = 0.2
        self.assertEqual(decay.scatter_fraction, 0.2)
        decay.scatter_fraction = 0.8
        self.assertEqual(decay.scatter_fraction, 0.8)

        decay.convolution_start = 12
        self.assertEqual(decay.convolution_start, 12)
        decay.convolution_start = 3
        self.assertEqual(decay.convolution_start, 3)

        decay.convolution_stop = 12
        self.assertEqual(decay.convolution_stop, 0)  # no data set
        decay.data = np.arange(20, dtype=np.float)
        self.assertEqual(decay.convolution_stop, 12)  # data set
        decay.convolution_stop = 3
        self.assertEqual(decay.convolution_stop, 3)

        decay.add_pile_up = True
        self.assertEqual(decay.add_pile_up, True)
        decay.add_pile_up = False
        self.assertEqual(decay.add_pile_up, False)

        decay.set_irf([1, 2, 3])
        self.assertListEqual(list(decay.get_irf()), [1, 2, 3])
        decay.set_irf([4, 5, 6])
        self.assertListEqual(list(decay.get_irf()), [4, 5, 6])

        decay.set_lifetime_spectrum([1, 2, 3, 4])
        self.assertListEqual(list(decay.get_lifetime_spectrum()), [1, 2, 3, 4])
        decay.set_lifetime_spectrum([4, 5, 6, 7])
        self.assertListEqual(list(decay.get_lifetime_spectrum()), [4, 5, 6, 7])

        decay.set_weights([1, 2, 3, 4])
        self.assertListEqual(list(decay.get_weights()), [1, 2, 3, 4])
        decay.set_weights([4, 5, 6, 7])
        self.assertListEqual(list(decay.get_weights()), [4, 5, 6, 7])

        decay.set_time_axis([1, 2, 3, 4])
        self.assertListEqual(list(decay.get_time_axis()), [1, 2, 3, 4])
        decay.set_time_axis([4, 5, 6, 7])
        self.assertListEqual(list(decay.get_time_axis()), [4, 5, 6, 7])

        decay.set_data([1, 2, 3, 45])
        self.assertListEqual(list(decay.get_data()), [1, 2, 3, 45])

        decay.irf_background = 892.1
        self.assertEqual(decay.irf_background, 892.1)

    def test_constructor_1(self):
        decay = fit2x.Decay()
        # default values
        self.assertEqual(decay.is_valid, False)
        self.assertEqual(decay.convolution_start, 0)
        self.assertEqual(decay.convolution_stop, 0)
        self.assertEqual(decay.add_pile_up, False)
        self.assertEqual(decay.use_amplitude_threshold, False)
        self.assertEqual(decay.excitation_period, 100.0)

    def test_constructor_2(self):
        decay = fit2x.Decay(
            decay_histogram=[1, 2, 3, 4, 56]
        )
        self.assertListEqual(list(decay.get_data()), [1, 2, 3, 4, 56])

        decay = fit2x.Decay(
            irf_histogram=[1., 2, 3, 4, 56]
        )
        self.assertListEqual(list(decay.get_irf()), [1, 2, 3, 4, 56])

        decay = fit2x.Decay(
            time_axis=[1., 2, 3, 4, 56]
        )
        self.assertListEqual(list(decay.get_time_axis()), [1, 2, 3, 4, 56])

        decay = fit2x.Decay(
            weights=[1., 2, 3, 4, 56]
        )
        self.assertListEqual(list(decay.get_weights()), [1, 2, 3, 4, 56])

        data = np.linspace(1, 22, 12)
        decay = fit2x.Decay(
            decay_histogram=data,
            convolution_range=(2, 32),
            add_pile_up=True,
            excitation_period=123.2
        )
        self.assertEqual(
            np.allclose(1. / np.sqrt(data), decay.get_weights()),
            True
        )
        self.assertEqual(len(decay.get_irf()), len(data))
        self.assertEqual(len(decay.get_time_axis()), len(data))
        self.assertEqual(decay.convolution_start, 2)
        self.assertEqual(decay.convolution_stop, min(len(data), 12))
        self.assertEqual(decay.add_pile_up, True)
        self.assertEqual(decay.excitation_period, 123.2)

    def test_convolve_lifetime_spectrum_variable_time_axis(self):
        time_axis = np.linspace(0, 25, 32)
        irf_position = 5.0
        irf_width = 1.0
        irf = scipy.stats.norm.pdf(time_axis, loc=irf_position, scale=irf_width)
        lifetime_spectrum = np.array([0.8, 1.1, 0.2, 4.0])
        model_decay = np.zeros_like(time_axis)

        fit2x.fconv_per_cs_time_axis(
            model_decay,
            time_axis=time_axis,
            lifetime_spectrum=lifetime_spectrum,
            instrument_response_function=irf
        )
        reference = np.array(
            [
                9.27891810e-07, 2.47480878e-05, 5.46047490e-04, 6.34298717e-03,
                3.99865961e-02, 1.41117009e-01, 2.92769678e-01, 3.83528048e-01,
                3.50102902e-01, 2.53687767e-01, 1.68661294e-01, 1.14128726e-01,
                8.12945833e-02, 6.06527275e-02, 4.67915552e-02, 3.69092130e-02,
                2.95268812e-02, 2.38266486e-02, 1.93277486e-02, 1.57274176e-02,
                1.28215151e-02, 1.04639958e-02, 8.54548307e-03, 6.98137666e-03,
                5.70483169e-03, 4.66231736e-03, 3.81060979e-03, 3.11463317e-03,
                2.54583913e-03, 2.08095096e-03, 1.70097037e-03, 1.39038160e-03
            ]
        )
        self.assertEqual(
            np.allclose(reference, model_decay, atol=1e-9),
            True
        )

    # def test_convolve_lifetime_spectrum_periodic(self):
    #     time_axis = np.linspace(0, 25, 25)
    #     irf_position = 6.0
    #     irf_width = 1.0
    #     irf = scipy.stats.norm.pdf(time_axis, loc=irf_position, scale=irf_width)
    #     lifetime_spectrum = np.array([0.2, 1.1, 0.8, 4.0])
    #     model_decay = np.zeros_like(time_axis)
    #     fit2x.Decay.convolve_lifetime_spectrum_periodic(
    #         model_decay,
    #         time_axis=time_axis,
    #         lifetime_spectrum=lifetime_spectrum,
    #         irf_histogram=irf,
    #         period=16.0
    #     )
    #     reference = np.array(
    #         [0.00560653, 0.00432208, 0.00342868, 0.00603476, 0.0454072,
    #          0.21058509, 0.4551221, 0.55543513, 0.48165047, 0.36885986,
    #          0.27946424, 0.21331303, 0.16359693, 0.12577494, 0.09681668,
    #          0.07457228, 0.05745678, 0.04427657, 0.03412254, 0.02629821,
    #          0.02026841, 0.01562132, 0.01203976, 0.00927939, 0.0071519])
    #     self.assertEqual(
    #         np.allclose(reference, model_decay), True
    #     )

    def test_shift(self):
        time_axis = np.linspace(0, 12, 25)
        irf_position = 6.0
        irf_width = 1.0
        irf = scipy.stats.norm.pdf(time_axis, loc=irf_position, scale=irf_width)
        # integer shift
        shifted_irf_ip = fit2x.Decay.shift_array(
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
        shifted_irf_in = fit2x.Decay.shift_array(
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
        shifted_irf_fp = fit2x.Decay.shift_array(
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

        shifted_irf_fn = fit2x.Decay.shift_array(
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
        shifted_irf_irp = fit2x.Decay.shift_array(
            input=irf,
            shift=26.0
        )
        self.assertEqual(
            np.allclose(shifted_irf_irp, shifted_irf_ip), True
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
        model_incl_irf = fit2x.Decay.add_curve(
            curve1=model_decay,
            curve2=irf,
            start=0,
            stop=-1,
            areal_fraction_curve2=0.9
        )
        ref = np.array(
            [0.34690072, 0.62173703, 1.01095837, 1.48560296, 1.97391665,
             2.37226495, 2.5796998, 2.53969438, 2.26563801, 1.83441674,
             1.352173, 0.91287367, 0.57136929, 0.33972356, 0.20072888,
             0.12605909, 0.08952876, 0.07268571, 0.06480473, 0.06056904,
             0.05769004, 0.05530237, 0.05311547, 0.05104114, 0.04905379,
             0.04714504, 0.04531079, 0.04354795, 0.04185369, 0.04022535,
             0.03866037, 0.03715627, 0.03571068, 0.03432134, 0.03298605,
             0.03170271, 0.0304693, 0.02928388, 0.02814458, 0.0270496,
             0.02599722, 0.02498578, 0.0240137, 0.02307943, 0.02218152,
             0.02131853, 0.02048913, 0.01969199, 0.01892586, 0.01818954,
             0.01748187, 0.01680173, 0.01614805, 0.0155198, 0.01491599,
             0.01433568, 0.01377794, 0.0132419, 0.01272672, 0.01223158,
             0.01175571, 0.01129834, 0.01085878, 0.01043631]
        )
        self.assertEqual(np.allclose(ref, model_incl_irf), True)

    def test_compute_decay(self):
        np.random.seed(0)
        time_axis = np.linspace(0, 10, 16)
        irf_position = 2.0
        irf_width = 0.5
        n_peak = 1000
        irf = scipy.stats.norm.pdf(
            time_axis,
            loc=irf_position,
            scale=irf_width
        )
        irf *= n_peak
        lifetime_spectrum = np.array([1.0, 4.])
        data_decay = np.zeros_like(time_axis)
        fit2x.fconv_per_cs_time_axis(
            data_decay,
            time_axis=time_axis,
            lifetime_spectrum=lifetime_spectrum,
            instrument_response_function=irf
        )
        data_decay = np.random.poisson(
            np.clip(data_decay, 1e-9, None)
        )
        data_weight = 1. / np.clip(data_decay, 1, 1e6)
        model = np.zeros_like(time_axis)
        irf += 0.0
        fit2x.Decay.compute_decay(
            model_function=model,
            data=data_decay,
            squared_weights=data_weight,
            time_axis=time_axis,
            irf_histogram=irf,
            lifetime_spectrum=lifetime_spectrum,
            scatter_fraction=0.1,
            excitation_period=5.,
            constant_offset=10,
            number_of_photons=-1,
            scale_model_to_data=False,
            use_amplitude_threshold=False
        )
        ref = np.array(
            [182.30544775, 150.93787478, 332.10434096, 752.30770162, 773.90059838,
             642.39549668, 543.50074262, 461.57573646, 392.25056234, 333.56811533,
             283.89449639, 241.84668574, 206.25398245, 176.12540958, 150.62212325,
             10.]
        )
        self.assertEqual(np.allclose(ref, model), True)

    def test_decay_class(self):
        time_axis, data = np.load('./data/reference/img_decay_histogram.npy').T
        data[0] = 0
        time_axis *= 4
        irf_position = 1.6
        irf_width = 0.08
        n_peak = 10000
        irf = scipy.stats.norm.pdf(
            time_axis,
            loc=irf_position,
            scale=irf_width
        )
        irf *= n_peak
        irf += 1e-6

        weights = 1. / (np.sqrt(data))
        weights[0] = 0
        decay_object = fit2x.Decay(
            decay_histogram=data,
            weights=weights,
            irf_histogram=irf,
            time_axis=time_axis,
            constant_offset=140.0,
            lifetime_spectrum=[1, 2]
        )
        m = decay_object.model
        wres = decay_object.weighted_residuals
        ref = np.array(
            [140., 140.00919308, 140.00118584, 140.00015296, 140.00001973,
             140.00000255, 140.00000033, 140.00000004, 140.00000001, 140.,
             140., 140., 140., ])

        # import pylab as p
        # p.plot(m[::64])
        # p.show()
        print(m[::64])
        self.assertEqual(
            np.allclose(m[::64], ref), True
        )
        wres_ref = np.array(
            [-0., 154.7475065, 54.30542351, 18.72900653, 1.20482741,
             -5.68219671, -7.8417569, -8.36660027, -6.25656282, -9.5,
             -9.90600991, -6.70820393, -7.17951632])
        self.assertEqual(
            np.allclose(wres[::64], wres_ref), True
        )

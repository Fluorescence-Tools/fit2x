from __future__ import division
import unittest

import scipy.stats
import fit2x
import numpy as np


class Tests(unittest.TestCase):

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
        reference = np.array([9.27891810e-07, 2.47480878e-05, 5.46047490e-04, 6.34298717e-03,
                              3.99865961e-02, 1.41117009e-01, 2.92769678e-01, 3.83528048e-01,
                              3.50102902e-01, 2.53687767e-01, 1.68661294e-01, 1.14128726e-01,
                              8.12945833e-02, 6.06527275e-02, 4.67915552e-02, 3.69092130e-02,
                              2.95268812e-02, 2.38266486e-02, 1.93277486e-02, 1.57274176e-02,
                              1.28215151e-02, 1.04639958e-02, 8.54548307e-03, 6.98137666e-03,
                              5.70483169e-03, 4.66231736e-03, 3.81060979e-03, 3.11463317e-03,
                              2.54583913e-03, 2.08095096e-03, 1.70097037e-03, 1.39038160e-03])
        self.assertEqual(
            np.allclose(reference, model_decay),
            True
        )

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
        time_axis = np.linspace(0, 10, 64)
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
        data_weight = 1. / np.clip(data_decay, 1, None)
        model = np.zeros_like(time_axis)
        irf += 0.0
        fit2x.Decay.compute_decay(
            model_function=model,
            data=data_decay,
            squared_weights=data_weight,
            time_axis=time_axis,
            instrument_response_function=irf,
            lifetime_spectrum=lifetime_spectrum,
            irf_background_counts=0.0,
            irf_shift_channels=-4.5,
            scatter_areal_fraction=0.1,
            excitation_period=5.,
            constant_background=10,
            total_area=-1,
            use_amplitude_threshold=False
        )
        ref = np.array(
            [200.14178357, 193.01541012, 186.6902567, 181.919255,
             180.2246366, 184.21582932, 197.621866, 224.62048753,
             268.24316442, 328.28031062, 399.87860574, 474.21554589,
             541.63228454, 595.78940302, 636.23678854, 667.56363864,
             695.68428619, 723.87049036, 751.10880506, 773.40438611,
             786.57771931, 788.51156644, 779.7515296, 762.70243287,
             740.34886113, 715.2820183, 689.29929699, 663.44266721,
             638.22297216, 613.84669358, 590.37261652, 567.79789123,
             546.09766628, 525.24077317, 505.19512305, 485.92931803,
             467.41305141, 449.61716818, 432.51364251, 416.07553804,
             400.27696623, 385.09304576, 370.4998633, 356.4744359,
             342.99467477, 330.03935048, 317.58805956, 305.62119233,
             294.11990206, 283.06607524, 272.44230308, 262.23185408,
             252.4186477, 242.98722903, 233.92274441, 225.2109181,
             216.83802974, 208.7908928, 201.05683373, 193.62367211,
             186.47970135, 179.61367034, 173.01476571, 166.67259476]
        )
        self.assertEqual(np.allclose(ref, model), True)

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

        decay.irf_shift_channels = 2.1
        self.assertAlmostEqual(decay.irf_shift_channels, 0.1)
        decay.irf_shift_channels = 2.2
        self.assertAlmostEqual(decay.irf_shift_channels, 0.2)

        # The number of photons is computed by scaling the model to the data
        # Here the data is zero. Hence, the number of photons is zero
        decay.number_of_photons = 11
        self.assertEqual(decay.number_of_photons, 0)

        # If the scaling is turned off we get the number of photons that
        # was specified by the attribute.
        decay.scale_model_to_data = False
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
        self.assertEqual(decay.convolution_stop, 1)  # no data set but minimal irf length is 1
        decay.data = np.arange(20, dtype=np.float)
        self.assertEqual(decay.convolution_stop, 12)  # data set
        decay.convolution_stop = 3
        self.assertEqual(decay.convolution_stop, 3)

        decay.use_pile_up_correction = True
        self.assertEqual(decay.use_pile_up_correction, True)
        decay.use_pile_up_correction = False
        self.assertEqual(decay.use_pile_up_correction, False)

        decay.set_irf([1, 2, 3])
        self.assertListEqual(list(decay.get_irf()), [1, 2, 3])
        decay.set_irf([4, 5, 6])
        self.assertListEqual(list(decay.get_irf()), [4, 5, 6])

        decay.set_lifetime_spectrum([1, 2, 3, 4])
        self.assertListEqual(list(decay.get_lifetime_spectrum()), [1, 2, 3, 4])
        decay.set_lifetime_spectrum([4, 5, 6, 7])
        self.assertListEqual(list(decay.get_lifetime_spectrum()), [4, 5, 6, 7])

        decay.set_data_weights([1, 2, 3, 4])
        self.assertListEqual(list(decay.get_data_weights()), [1, 2, 3, 4])
        decay.set_data_weights([4, 5, 6, 7])
        self.assertListEqual(list(decay.get_data_weights()), [4, 5, 6, 7])

        decay.set_time_axis([1, 2, 3, 4])
        self.assertListEqual(list(decay.get_time_axis()), [1, 2, 3, 4])
        decay.set_time_axis([4, 5, 6, 7])
        self.assertListEqual(list(decay.get_time_axis()), [4, 5, 6, 7])

        decay.set_data([1, 2, 3, 45])
        self.assertListEqual(list(decay.get_data()), [1, 2, 3, 45])

        decay.irf_background = 892.1
        self.assertEqual(decay.irf_background, 892.1)

    def test_parameter(self):
        decay = fit2x.Decay()
        decay.lifetime_spectrum = [1., 4.]
        ref = {'irf_background_counts': 0.0,
               'irf_shift_channels': 0.0,
               'scatter_fraction': 0.0,
               'constant_offset': 0.0,
               'lifetime_spectrum': np.array([1., 4.], dtype=np.float64),
               'number_of_photons': 0.0,
               'acquisition_time': np.finfo(float).max,
               'instrument_dead_time': np.finfo(float).eps,
               'convolution_range': (0, -1),
               'excitation_period': np.finfo(float).max,
               'scale_model_to_data': True,
               'score_range': (0, -1),
               'use_corrected_irf_as_scatter': False,
               'amplitude_threshold': np.finfo(float).eps,
               'use_amplitude_threshold': False,
               'use_pile_up_correction': False,
               'convolution_method': 0,
               'use_linearization': False}
        a = decay.parameter
        # test lifetime spectrum separately as == is not well defined for np.array
        self.assertEqual(
            np.alltrue(a.pop('lifetime_spectrum') == ref.pop('lifetime_spectrum')),
            True
        )
        self.assertDictEqual(a, ref)
        ref['irf_background_counts'] = 11
        decay.set(**ref)
        self.assertEqual(decay.parameter['irf_background_counts'], ref['irf_background_counts'])

    def test_constructor_1(self):
        decay = fit2x.Decay()
        # default values
        self.assertEqual(decay.is_valid, False)
        self.assertEqual(decay.convolution_start, 0)
        self.assertEqual(decay.convolution_stop, 0)
        self.assertEqual(decay.use_pile_up_correction, False)
        self.assertEqual(decay.use_amplitude_threshold, False)
        self.assertEqual(decay.excitation_period, np.finfo(float).max)

    def test_constructor_2(self):
        decay = fit2x.Decay(
            data=[1, 2, 3, 4, 56]
        )
        self.assertListEqual(list(decay.data), [1, 2, 3, 4, 56])

        decay = fit2x.Decay(
            irf_histogram=[1., 2, 3, 4, 56]
        )
        self.assertListEqual(list(decay.irf), [1, 2, 3, 4, 56])

        decay = fit2x.Decay(
            time_axis=[1., 2, 3, 4, 56]
        )
        self.assertListEqual(list(decay.time_axis), [1, 2, 3, 4, 56])

        decay = fit2x.Decay(
            data_weights=[1., 2, 3, 4, 56]
        )
        self.assertListEqual(list(decay.data_weights), [1, 2, 3, 4, 56])

        data = np.linspace(1, 22, 12)
        decay = fit2x.Decay(
            data=data,
            convolution_range=(2, 32),
            use_pile_up_correction=True,
            excitation_period=123.2
        )
        self.assertEqual(
            np.allclose(1. / np.sqrt(data), decay.data_weights),
            True
        )
        self.assertEqual(len(decay.get_irf()), len(data))
        self.assertEqual(len(decay.get_time_axis()), len(data))
        self.assertEqual(decay.convolution_start, 2)
        self.assertEqual(decay.convolution_stop, min(len(data), 12))
        self.assertEqual(decay.use_pile_up_correction, True)
        self.assertEqual(decay.excitation_period, 123.2)

    def test_convolve_lifetime_spectrum_variable_time_axis(self):
        period = 25
        time_axis = np.linspace(0, period, 16)
        irf_position = 5.0
        irf_width = 1.0
        irf = scipy.stats.norm.pdf(time_axis, loc=irf_position, scale=irf_width)
        irf /= np.max(irf)
        lifetime_spectrum = np.array([0.8, 1.1, 0.2, 4.0])
        model_decay = np.zeros_like(time_axis)
        irf[irf < 0.0001] = 0.0

        fit2x.fconv_per_cs_time_axis(
            model_decay,
            convolution_stop=len(irf),
            convolution_start=0,
            time_axis=time_axis,
            lifetime_spectrum=lifetime_spectrum,
            instrument_response_function=irf,
            period=period
        )
        print(model_decay)
        reference = np.array(
            [0.00349525, 0.00552579, 0.21129489, 0.96300658, 0.77383893, 0.36838528,
             0.18074355, 0.10503629, 0.06614159, 0.04292135, 0.02814564, 0.01852181,
             0.01220309, 0.00804318, 0.00530204, 0.00349525]
        )
        # import pylab as plt
        # plt.semilogy(irf)
        # plt.semilogy(reference)
        # plt.semilogy(model_decay)
        # plt.show()
        self.assertEqual(
            np.allclose(reference, model_decay, atol=1e-9),
            True
        )

    def test_convolve_lifetime_spectrum_periodic(self):
        period = 25
        time_axis = np.linspace(0, period, 64)
        irf_position = 6.0
        irf_width = 0.5
        irf = scipy.stats.norm.pdf(time_axis, loc=irf_position, scale=irf_width)
        irf /= np.max(irf)
        irf[irf < 0.0001] = 0.0
        lifetime_spectrum = np.array([0.2, 1.1 , 0.8, 4.1])
        model_decay = np.zeros_like(time_axis)
        fit2x.fconv_per_cs_time_axis(
            model_decay,
            time_axis=time_axis,
            lifetime_spectrum=lifetime_spectrum,
            instrument_response_function=irf,
            convolution_start=0,
            convolution_stop=len(irf),
            period=period
        )
        reference = np.array([0.00987956, 0.78140209, 0.19914572, 0.04219551])
        print(model_decay[::16])
        # import pylab as plt
        # plt.semilogy(irf, label='irf')
        # plt.semilogy(reference, label='ref')
        # plt.semilogy(model_decay, label='model')
        # plt.legend()
        # plt.show()
        self.assertEqual(
            np.allclose(reference, model_decay[::16]), True
        )

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
            instrument_response_function=irf,
            convolution_stop=len(irf),
            convolution_start=1
        )
        areal_fraction_2 = 0.9
        ref = model_decay / np.sum(model_decay) * (1-areal_fraction_2) + irf / np.sum(irf) * areal_fraction_2
        ref *= np.sum(model_decay)
        model_incl_irf = fit2x.Decay.add_curve(
            curve1=model_decay,
            curve2=irf,
            start=0,
            stop=-1,
            areal_fraction_curve2=areal_fraction_2
        )
        self.assertEqual(np.allclose(ref, model_incl_irf), True)

    def test_compute_decay(self):
        np.random.seed(0)
        period = 13.6
        time_axis = np.linspace(0, period, 64)
        irf_position = 2.0
        irf_width = 0.5
        n_peak = 1000
        irf = scipy.stats.norm.pdf(
            time_axis,
            loc=irf_position,
            scale=irf_width
        )
        irf *= n_peak
        irf[irf < 0.01] = 0.0
        lifetime_spectrum = np.array([1.0, 4.])
        data = np.zeros_like(time_axis)
        fit2x.fconv_per_cs_time_axis(
            data,
            time_axis=time_axis,
            lifetime_spectrum=lifetime_spectrum,
            instrument_response_function=irf,
            convolution_start=1,
            convolution_stop=len(irf),
            period=period
        )
        data = np.random.poisson(
            np.clip(data, 1e-9, None)
        )
        data_weight = 1. / np.clip(data, 1, 1e6)
        model = np.zeros_like(time_axis)
        irf += 0.0
        fit2x.Decay.compute_decay(
            model_function=model,
            data=data,
            squared_data_weights=data_weight,
            time_axis=time_axis,
            irf_histogram=irf,
            lifetime_spectrum=lifetime_spectrum,
            scatter=irf,
            scatter_fraction=0.1,
            excitation_period=period,
            constant_offset=10,
            number_of_photons=-1,
            scale_model_to_data=True,
            use_amplitude_threshold=False,
            linearization=np.ones_like(model),
            convolution_start=0,
            convolution_stop=len(irf),
        )

        # print(model[::16])
        # import pylab as plt
        # plt.semilogy(time_axis, irf)
        # plt.semilogy(time_axis, model)
        # plt.show()

        ref = np.array([ 61.73859508, 665.13216763, 285.10189626, 126.00687973])
        self.assertEqual(np.allclose(ref, model[::16]), True)

    def test_decay_class(self):
        time_axis, data = np.load('./data/reference/img_decay_histogram.npy').T
        data[0] = 0
        time_axis *= 4
        irf_position = 12.6
        irf_width = 0.1
        n_peak = 10000
        irf = scipy.stats.norm.pdf(
            time_axis,
            loc=irf_position,
            scale=irf_width
        )
        irf *= n_peak
        irf[irf < 1] = 0.0
        #irf += 1e-6

        weights = 1. / (np.sqrt(np.clip(data, 1, 1e12)))
        weights[0] = 0
        period = np.max(time_axis)
        print("time_axis:", time_axis.shape)
        print("data:", data.shape)
        print("irf:", irf.shape)
        print("weights:", weights.shape)

        lt = [1, 0.2, 0.1, 4]
        decay_object = fit2x.Decay(
            data=data,
            data_weights=weights,
            irf_histogram=irf,
            time_axis=time_axis,
            constant_offset=0.0,
            lifetime_spectrum=lt,
            scale_model_to_data=False,
            convolution_range=(0, len(irf)),
            excitation_period=period
        )
        m = decay_object.model
        wres = decay_object.weighted_residuals
        ref = np.array([6.53340861e-07, 4.13128530e-03, 6.87408535e-05, 1.14378588e-06])

        # import pylab as plt
        # plt.semilogy(decay_object.time_axis, decay_object.model)
        # plt.semilogy(decay_object.time_axis, decay_object.irf)
        # plt.show()

        self.assertEqual(np.allclose(m[::256], ref), True)
        wres_ref = np.array([-0.0 , 12.44956777,  9.11042603,  8.77496426])
        self.assertEqual(np.allclose(wres[::256], wres_ref), True)

        parameter = decay_object.parameter
        parameter['score_range'] = (0, 64)
        parameter['lifetime_spectrum'] = [1, 0.2, 0.1, 4]
        score = decay_object.compute_score(
            data=decay_object.data,
            time_axis=decay_object.time_axis,
            irf_histogram=decay_object.irf,
            data_weights=decay_object.data_weights,
            **parameter
        )
        decay_object.lifetime_spectrum = [1, 0.2, 0.1, 4]
        score_decay = decay_object.get_score(0, 64, 'poisson')
        self.assertAlmostEqual(score_decay, score)

    # def test_convolution_methods(self):
    #     from matplotlib import pylab as plt
    #     import scipy.stats
    #     import fit2x
    #
    #     decay = fit2x.Decay()
    #     decay.lifetime_spectrum = [1., 0.5, 0.1, 4]
    #     period = 25
    #     time_axis = np.linspace(0, period, 16)
    #     irf_position = 5.0
    #     irf_width = 0.3
    #     irf = scipy.stats.norm.pdf(time_axis, loc=irf_position, scale=irf_width)
    #     decay.irf = irf
    #     decay.scale_model_to_data = False
    #
    #
    #     plt.semilogy(irf)
    #
    #     decay.convolution_method = 0
    #     plt.semilogy(decay.model, label=decay.convolution_method)
    #
    #     decay.convolution_method = 1
    #     plt.semilogy(decay.model, label=decay.convolution_method)
    #
    #     decay.convolution_method = 2
    #     plt.semilogy(decay.model, label=decay.convolution_method)
    #
    #     decay.convolution_method = 3
    #     plt.semilogy(decay.model, label=decay.convolution_method)
    #
    #     decay.convolution_method = 4
    #     plt.semilogy(decay.model, label=decay.convolution_method)
    #
    #     decay.convolution_method = 5
    #     plt.semilogy(decay.model, label=decay.convolution_method)
    #
    #     plt.legend()
    #     plt.show()
    #
    #     a = decay.parameter
    #     # test lifetime spectrum separately as == is not well defined for np.array
    #     self.assertEqual(
    #         np.alltrue(a.pop('lifetime_spectrum') == ref.pop('lifetime_spectrum')),
    #         True
    #     )
    #     self.assertDictEqual(a, ref)
    #     ref['irf_background_counts'] = 11
    #     decay.set(**ref)
    #     self.assertEqual(decay.parameter['irf_background_counts'], ref['irf_background_counts'])
    #
    #

from __future__ import division
import unittest

import scipy.stats
import fit2x
import numpy as np


class Tests(unittest.TestCase):

    # # OK - reference needs update
    # def test_convolve_lifetime_spectrum_variable_time_axis(self):
    #     time_axis = np.linspace(0, 25, 32)
    #     irf_position = 5.0
    #     irf_width = 1.0
    #     irf = scipy.stats.norm.pdf(time_axis, loc=irf_position, scale=irf_width)
    #     lifetime_spectrum = np.array([0.8, 1.1, 0.2, 4.0])
    #     model_decay = np.zeros_like(time_axis)
    #     fit2x.fconv_per_cs_time_axis(
    #         model_decay,
    #         time_axis=time_axis,
    #         lifetime_spectrum=lifetime_spectrum,
    #         instrument_response_function=irf
    #     )
    #     reference = np.array([9.27891810e-07, 2.47480878e-05, 5.46047490e-04, 6.34298717e-03,
    #                           3.99865961e-02, 1.41117009e-01, 2.92769678e-01, 3.83528048e-01,
    #                           3.50102902e-01, 2.53687767e-01, 1.68661294e-01, 1.14128726e-01,
    #                           8.12945833e-02, 6.06527275e-02, 4.67915552e-02, 3.69092130e-02,
    #                           2.95268812e-02, 2.38266486e-02, 1.93277486e-02, 1.57274176e-02,
    #                           1.28215151e-02, 1.04639958e-02, 8.54548307e-03, 6.98137666e-03,
    #                           5.70483169e-03, 4.66231736e-03, 3.81060979e-03, 3.11463317e-03,
    #                           2.54583913e-03, 2.08095096e-03, 1.70097037e-03, 1.39038160e-03])
    #     self.assertEqual(
    #         np.allclose(reference, model_decay),
    #         True
    #     )

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

    def test_getter_setter(self):
        decay = fit2x.Decay(
            scale_model_to_data=True
        )

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
        # If there is not data convolution starts at zero
        self.assertEqual(decay.convolution_start, 0)
        decay.data = np.ones(200)
        self.assertEqual(decay.convolution_start, 12)

        decay.convolution_stop = 12
        decay.data = []
        self.assertEqual(decay.convolution_stop, 0)
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
        ref = {
               'acquisition_time': 1000000000.0,
               'use_amplitude_threshold': False,
               'abs_lifetime_spectrum': False,
               'amplitude_threshold': 2.220446049250313e-16,
               'convolution_range': (0, 0),
               'use_corrected_irf_as_scatter': True,
               'scatter_fraction': 0.0,
               'convolution_method': 0,
               'excitation_period': 100.0,
               'irf_shift_channels': 0.0,
               'irf_background_counts': 0.0,
               'constant_offset': 0.0,
               'pile_up_model': 'coates',
               'instrument_dead_time': 120.0,
               'use_pile_up_correction': False,
               'scale_model_to_data': False,
               'number_of_photons': -1.0,
               'data': np.array([], dtype=np.float64),
               'linearization_table': np.array([], dtype=np.float64),
               'data_weights': np.array([], dtype=np.float64),
               'time_axis': np.array([], dtype=np.float64),
               'irf_histogram': np.array([], dtype=np.float64),
               'lifetime_spectrum': np.array([1., 4.]),
               'use_linearization': False,
               'score_range': (0, -1),
               'score_type': 'poisson'}
        a = decay.parameter
        # test lifetime spectrum separately as == is not well defined for np.array
        for n in [
            'data', 'linearization_table',
            'data_weights', 'time_axis',
            'irf_histogram', 'lifetime_spectrum']:
            self.assertEqual(
                np.alltrue(a.pop(n) == ref.pop(n)),
                True
            )
        self.assertDictEqual(a, ref)

    def test_constructor_1(self):
        decay = fit2x.Decay()
        # default values
        self.assertEqual(decay.is_valid, False)
        self.assertEqual(decay.convolution_start, 0)
        self.assertEqual(decay.convolution_stop, 0)
        self.assertEqual(decay.use_pile_up_correction, False)
        self.assertEqual(decay.use_amplitude_threshold, False)
        self.assertEqual(decay.excitation_period, 100.0)

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


    def test_compute_decay(self):
        # import numpy as np
        # import fit2x
        # import scipy.stats
        np.random.seed(0)
        period = 13.6
        time_axis = np.linspace(0, period, 16)
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
        lifetime_spectrum = np.array([1.0, 4., 2., 0.1])
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
        data_weight = 1./np.clip(np.sqrt(data), 1, 1e9)
        irf += 0.0
        decay = fit2x.Decay(
            data=data,
            data_weights=data_weight,
            time_axis=time_axis,
            irf_histogram=irf,
            lifetime_spectrum=lifetime_spectrum,
            scatter_fraction=0.1,
            excitation_period=period,
            constant_offset=10,
            number_of_photons=-1,
            scale_model_to_data=False,
            use_amplitude_threshold=False,
            convolution_range=(0, -1),
        )
        model = decay.model

        # print(model[::16])
        # import pylab as plt
        # plt.semilogy(time_axis, irf)
        # plt.semilogy(time_axis, model)
        # plt.show()

        ref = np.array([  56.83308935,  159.6799803 , 1229.85116673,  937.6613214 ,
                          572.69332145,  455.9380128 ,  365.4952043 ,  293.39596583,
                          235.91942871,  190.09991116,  153.57321186,  124.45462149,
                          101.24167532,   82.73662878,   67.9846561 ,   59.87080939])
        self.assertEqual(np.allclose(ref, model), True)


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
        self.assertEqual(decay.convolution_stop, 0)  # no data set
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
            data_weights=[1., 2, 3, 4, 56]
        )
        self.assertListEqual(list(decay.get_data_weights()), [1, 2, 3, 4, 56])

        data = np.linspace(1, 22, 12)
        decay = fit2x.Decay(
            data=data,
            convolution_range=(2, 32),
            use_pile_up_correction=True,
            excitation_period=123.2
        )
        self.assertEqual(
            np.allclose(1. / np.sqrt(data), decay.get_data_weights()),
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
        data_decay = np.zeros_like(time_axis)
        fit2x.fconv_per_cs_time_axis(
            data_decay,
            time_axis=time_axis,
            lifetime_spectrum=lifetime_spectrum,
            instrument_response_function=irf,
            convolution_start=1,
            convolution_stop=len(irf),
            period=period
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
        print(model[::16])
        # import pylab as plt
        # plt.semilogy(time_axis, irf)
        # plt.semilogy(time_axis, model)
        # plt.show()
        ref = np.array([59.72502663, 639.31947171, 274.26263987, 121.43610676])
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

        decay_object = fit2x.Decay(
            data=data,
            data_weights=weights,
            irf_histogram=irf,
            time_axis=time_axis,
            constant_offset=0.0,
            lifetime_spectrum=[1, 2],
            scale_model_to_data=False,
            convolution_range=(0, len(irf)),
            excitation_period=period
        )
        m = decay_object.model
        wres = decay_object.weighted_residuals
        ref = np.array([1.20803660e-10, 1.55827821e-11, 2.01006408e-12, 3.86467220e-05,
                        4.83026185e-03, 6.23068189e-04, 8.03712057e-05, 1.03672934e-05,
                        1.33730446e-06, 1.72502422e-07, 2.22515415e-08, 2.87028490e-09,
                        3.70245605e-10])

        # import pylab as plt
        # plt.semilogy(decay_object.time_axis, decay_object.model)
        # plt.semilogy(decay_object.time_axis, decay_object.irf)
        # plt.show()

        self.assertEqual(
            np.allclose(m[::64], ref), True
        )
        wres_ref = np.array([ -0.        , 155.64703659,  56.77147171,  24.45403694,
                              12.44951162,   9.32731225,   8.54399434,   8.36659903,
                              9.11043343,   7.99999998,   7.87400787,   8.94427191,
                              8.77496439])
        self.assertEqual(
            np.allclose(wres[::64], wres_ref), True
        )

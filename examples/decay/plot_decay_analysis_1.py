"""
===============================
Fluorescence decay analysis - 1
===============================
In many cases fluorescence decays can be described by linear combinations of
exponential decays. The ``Decay`` class of ``tttrlib`` computes models for fluorescence
decays, f(t), fluorescence intensity decay histograms, h(t), and to scores models
against experimental decay histograms. The ``Decay`` class handles typical artifacts
encountered in fluorescence decays such as scattered light, constant backgrounds,
convolutions with the instrument response function :cite:`oconnor_time-correlated_2012`,
pile-up artifacts :cite:`coates_correction_1968`, differential non-linearity of the
micro time channels :cite:`becker_advanced_2005` and more.

Below the basic usage of the ``Decay`` class is outlined with a few application
examples. These examples can be used as starting point for custom analysis pipelines,
e.g. for burst-wise single-molecule analysis, pixel-wise FLIM analysis or analysis
over larger ensembles of molecules or pixels.

Decay histograms
----------------

.. code-block::python

    tttr_object = tttrlib.TTTR('../../test/data/BH/BH_SPC132.spc', 'SPC-130')
    micro_time_coarsening = 8
    counts, bins = fit2x.Decay.compute_microtime_histogram(
        tttr_object,
        micro_time_coarsening=micro_time_coarsening,
    )



Decay class
-----------
Fluorescence decays can be either computed using the static method provided by the
``Decay`` class or
Create an instance the ``Decay`` class

.. code-block::python

    decay = fit2x.Decay(
        data=data_decay.astype(np.float64),
        instrument_response_function=irf.astype(np.float64),
        time_axis=time_axis,
        period=data.header.macro_time_resolution
    )


Single-molecule

"""
import matplotlib.pylab as plt
import scipy.optimize
import scipy.stats
import numpy as np
import tttrlib
import fit2x


def objective_function_chi2(
        x: np.ndarray,
        decay: fit2x.Decay,
        x_min: int = 20,
        x_max: int = 150
):
    scatter, background, time_shift, irf_background = x[0:4]
    lifetime_spectrum = np.abs(x[4:])
    decay.set_lifetime_spectrum(lifetime_spectrum)
    decay.irf_background_counts = irf_background
    decay.scatter_fraction = scatter
    decay.constant_offset = background
    decay.irf_shift = time_shift
    # wres = decay.get_weighted_residuals()
    # return np.sum(wres[x_min:x_max]**2)
    return decay.get_score(x_min, x_max)


def objective_function_mle(
        x: np.ndarray,
        decay: fit2x.Decay,
        x_min: int = 20,
        x_max: int = 500
):
    scatter, background, time_shift, irf_background = x[0:4]
    lifetime_spectrum = np.abs(x[4:])
    decay.set_lifetime_spectrum(lifetime_spectrum)
    decay.irf_background_counts = irf_background
    decay.scatter_fraction = scatter
    decay.constant_offset = background
    decay.irf_shift = time_shift
    chi2_mle = decay.get_score(x_min, x_max, score_type="poisson")
    return chi2_mle


#######################################
#           Prepare data
#######################################

# load file
spc132_filename = '../../test/data/bh_spc132_sm_dna/m000.spc'
data = tttrlib.TTTR(spc132_filename, 'SPC-130')
ch0_indeces = data.get_selection_by_channel([0, 8])
data_ch0 = data[ch0_indeces]

n_bins = 128
# selection from tttr object
cr_selection = data_ch0.get_selection_by_count_rate(
    time_window=6.0e-3, n_ph_max=7
)
low_count_selection = data_ch0[cr_selection]
# create histogram for IRF
irf, _ = np.histogram(low_count_selection.micro_times, bins=n_bins)

# Select high count regions
# equivalent selection using selection function
cr_selection = tttrlib.selection_by_count_rate(
    data_ch0.macro_times,
    0.100e-3, n_ph_max=5,
    macro_time_calibration=data.header.macro_time_resolution,
    invert=True
)
high_count_selection = data_ch0[cr_selection]
data_decay, _ = np.histogram(high_count_selection.micro_times, bins=n_bins)
time_axis = np.arange(0, n_bins) * data.header.micro_time_resolution * 4096 * 1e9 / n_bins  # in [ns]

###################################
#     Make decay object
###################################
decay = fit2x.Decay(
    data=data_decay.astype(np.float64),
    irf_histogram=irf.astype(np.float64),
    time_axis=time_axis,
    scale_model_to_data=True,
    excitation_period=max(time_axis),
    lifetime_spectrum=[1., 1.2, 1., 3.5],
    abs_lifetime_spectrum=True,
    convolution_method=fit2x.ConvFastPeriodicTime
)
decay.update_model()

# A minimum number of photons should be in each channel
# as no MLE is used and Gaussian distributed errors are assumed
sel = np.where(data_decay > 1)[0]
x_min = 2 #int(min(sel))

# The BH card are not linear at the end of the TAC. Thus
# fit not the full range
x_max = 110 #max(sel)

# Set some initial values for the fit
scatter = [0.05]
background = [0.01]
time_shift = [0]
irf_background = [5]
lifetime_spectrum = [0.5, 0.5, 0.5, 2.5]

###################################
#     Optimize                    #
###################################
x0 = np.hstack(
    [
        scatter,
        background,
        time_shift,
        irf_background,
        lifetime_spectrum
    ]
)


x0 = np.hstack([
    scatter, background, 
    time_shift, irf_background, 
    lifetime_spectrum])
fit = scipy.optimize.minimize(
    objective_function_chi2, x0,
    args=(decay, x_min, x_max)
)
fit_map = fit.x

###################################
#     Plot                        #
###################################
fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax[1].semilogy(decay.time_axis, decay.corrected_irf, label="IRF")
ax[1].semilogy(decay.time_axis, decay.data, label="Data")
ax[1].semilogy(decay.time_axis[x_min:x_max], decay.model[x_min:x_max], label="Model")
ax[1].legend()
ax[0].plot(decay.time_axis[x_min:x_max], 
    decay.weighted_residuals[x_min:x_max], 
    label='w.res.', color='green')
plt.show()


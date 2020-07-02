
// File: index.xml

// File: classap_1_1ap__error.xml


%feature("docstring") ap::ap_error "
";

%feature("docstring") ap::ap_error::make_assertion "
";

// File: classbfgs.xml


%feature("docstring") bfgs "
";

%feature("docstring") bfgs::bfgs "
";

%feature("docstring") bfgs::bfgs "
";

%feature("docstring") bfgs::~bfgs "
";

%feature("docstring") bfgs::setN "
";

%feature("docstring") bfgs::seteps "
";

%feature("docstring") bfgs::seteps "
";

%feature("docstring") bfgs::fix "
";

%feature("docstring") bfgs::free "
";

%feature("docstring") bfgs::minimize "
";

// File: classap_1_1complex.xml


%feature("docstring") ap::complex "
";

%feature("docstring") ap::complex::complex "
";

%feature("docstring") ap::complex::complex "
";

%feature("docstring") ap::complex::complex "
";

%feature("docstring") ap::complex::complex "
";

// File: classap_1_1const__raw__vector.xml


%feature("docstring") ap::const_raw_vector "
";

%feature("docstring") ap::const_raw_vector::const_raw_vector "
";

%feature("docstring") ap::const_raw_vector::GetData "
";

%feature("docstring") ap::const_raw_vector::GetLength "
";

%feature("docstring") ap::const_raw_vector::GetStep "
";

// File: class_decay.xml


%feature("docstring") Decay "
";

%feature("docstring") Decay::compute_decay "

Computes a fluorescence decay for a lifetime spectrum.  

The lifetime spectrum is a interleaved array of amplitudes and fluorescence
lifetimes. The fluorescence decay for the lifetime spectrum is computed. The
computed decay is convolved with a instrument response function (IRF). Before
convolution the IRF is corrected for a constant offset. The convolution
considers periodic excitation. The IRF is shifted by a specified value of micro
time channels. After convolution a constant fraction of scattered light is added
to the computed decay and the decay is scaled to the number of photons in an
experimental fluorescence decay. Finally, a constant background is added.
Optionally, pile-up is added to the computed fluorescence decay. Before the
fluorescence decay is computed an amplitude threshold can be applied to the
fluorescence lifetime spectrum to discriminate fluorescence lifetimes with small
amplitudes.  

Parameters
----------
* `model_function[out]` :  
    The output array that will contain the computed fluorescence decay  
* `n_model_function[out]` :  
    The number of points in the model function  
* `data[in]` :  
    The data for which the model function is computed. The data will be used to
    scaled the computed decay.  
* `n_data[in]` :  
    The number of data points  
* `squared_weights[in]` :  
    The squared weights of the data points. The data weights are used to scale
    the model function to the data (usually the data weights is Poissonian)  
* `n_weights[in]` :  
    The number of weights  
* `time_axis[in]` :  
    The time axiss for which the fluorescence decay is computed for  
* `n_time_axis[in]` :  
    The number of points in the time axis  
* `instrument_response_function[in]` :  
    The instrument response function used for convolution  
* `n_instrument_response_function[in]` :  
    The number of points in the instrument response function  
* `lifetime_spectrum[in]` :  
    The lifetime spectrum. A lifetime specturm is an interleaved array of
    amplitudes and fluorescence lifetimes, .e.g, (amplitude1, lifetime1,
    amplitude2, lifetime2, ...)  
* `n_lifetime_spectrum[in]` :  
    The number of points in the fluorescence lifetime spectrum  
* `start[in]` :  
    The start of the convolution  
* `stop[in]` :  
    The stop of the convoution. The decay will be computed in the range [start,
    stop]  
* `irf_background_counts[in]` :  
    The background counts of the instrument response function. This number will
    be subtracted from the IRF before convolution.  
* `irf_shift_channels[in]` :  
    The number of micro time channels the IRF will be shifted before the
    fluorescence lifetimes are convoluted with the IRF.  
* `scatter_areal_fraction[in]` :  
    The fraction (integrated fraction), i.e., the area the scattered light will
    have in the computed decay.  
* `excitation_period[in]` :  
    The excitation period (in units of the fluorescence lifetime, usually
    nanoseconds) that was used to excite the sample  
* `constant_background` :  
    A constant offset that is added to the fluorescence decay.  
* `total_area[in]` :  
    the area to which the model fluorescence decay is scaled to. If this value
    is negative (default value is -1) the model fluorescence decay is scaled to
    the data in the range defined by the parameters start, stop  
* `use_amplitude_threshold[in]` :  
    if set to true (default is false) a threshold will be applied. If this
    parameter is set to true, fluorescence lifetimes with a amplitude that is
    smaller than a threshold will be not considered.  
* `amplitude_threshold` :  
    The threshold that is used to discriminate fluorescence lifetimes that are
    smaller.  
* `pile_up` :  
    if set to true (default is false) pile up will be added to the model
    function.  
* `instrument_dead_time` :  
    the dead time of the instrument (used for pile up (in units of the lifetime,
    usually nano seconds)  
* `acquisition_time` :  
    the total time the acquisition of the decay in seconds.  
* `add_corrected_irf` :  
    if set to true (default is false) the background corrected irf will be added
    as scatter to the decay. If this is false the irf prior to a background and
    shift corrected irf is added as a scatter fraction.  
* `scale_model_to_area` :  
    if set to true (default is true) the model is either scaled to the area s
    provided by total_area (if total_area is larger then zero) or to the total
    number of experimental counts.  
";

%feature("docstring") Decay::shift_array "

Shift an input array by a floating number.  

Parameters
----------
* `input[in]` :  
    the input array  
* `n_input[in]` :  
    length of the input array  
* `output` :  
    output array  
* `n_output` :  
    length of the output array  
* `shift[in]` :  
    the shift of the output  
";

%feature("docstring") Decay::add_curve "

Computes the sum of two arrays considering their respective areal fraction.  

A weighted sum of two curves is computed. The weighted sum is computed by the
area of the first curve and the areal fraction of the second curve. The area of
the computed curve equals to the area of the first input curve while the areal
fraction of the second input curve will be equal to the specified value in the
range specified by the input parameters.  

Parameters
----------
* `output` :  
    the computed output curve (array)  
* `n_output` :  
    the number of points in the output  
* `curve1[in]` :  
    the first input curve / array  
* `n_curve1[in]` :  
    number of points in the first array  
* `curve2[in]` :  
    second curve / array  
* `n_curve2[in]` :  
    number of points in the second array  
* `areal_fraction_curve2[in]` :  
    areal fraction of the second curve in the output array  
* `start[in]` :  
    start index used for the area calculation  
* `stop[in]` :  
    stop index used for the area calculation  
";

%feature("docstring") Decay::is_valid "
";

%feature("docstring") Decay::set_data "
";

%feature("docstring") Decay::get_data "
";

%feature("docstring") Decay::set_use_amplitude_threshold "
";

%feature("docstring") Decay::get_use_amplitude_threshold "
";

%feature("docstring") Decay::set_amplitude_threshold "
";

%feature("docstring") Decay::get_amplitude_threshold "
";

%feature("docstring") Decay::set_total_area "
";

%feature("docstring") Decay::get_total_area "
";

%feature("docstring") Decay::set_period "
";

%feature("docstring") Decay::set_score_range "
";

%feature("docstring") Decay::get_period "
";

%feature("docstring") Decay::set_irf_shift_channels "
";

%feature("docstring") Decay::get_irf_shift_channels "
";

%feature("docstring") Decay::set_areal_scatter_fraction "
";

%feature("docstring") Decay::get_areal_scatter_fraction "
";

%feature("docstring") Decay::set_constant_background "
";

%feature("docstring") Decay::get_constant_background "
";

%feature("docstring") Decay::set_convolution_start "
";

%feature("docstring") Decay::get_convolution_start "
";

%feature("docstring") Decay::set_convolution_stop "
";

%feature("docstring") Decay::get_convolution_stop "
";

%feature("docstring") Decay::set_correct_pile_up "
";

%feature("docstring") Decay::get_correct_pile_up "
";

%feature("docstring") Decay::set_irf "
";

%feature("docstring") Decay::get_irf "
";

%feature("docstring") Decay::get_model "
";

%feature("docstring") Decay::set_lifetime_spectrum "
";

%feature("docstring") Decay::get_lifetime_spectrum "
";

%feature("docstring") Decay::set_weights "
";

%feature("docstring") Decay::get_weights "
";

%feature("docstring") Decay::set_time_axis "
";

%feature("docstring") Decay::get_time_axis "
";

%feature("docstring") Decay::set_irf_background_counts "
";

%feature("docstring") Decay::get_irf_background_counts "
";

%feature("docstring") Decay::Decay "

Parameters
----------
* `tttr_data` :  
    pointer to TTTR object that is used to construct a decay histogram  
* `micro_time_coarsening` :  
    an (optional) integer by which the micro times are divided to coarsen the
    time axis (default is 1)  
* `decay_data` :  
    the data to which the decay is fitted  
* `time_axis` :  
    the time axis that belongs to the data  
* `dt` :  
    the spacing between the points in the time axis. This optional parameter is
    used to compute a time axis if not time axis was provided by the parameter
    time_axis  
* `weights` :  
    the weights of the data points. If the weights are not provided (nullptr /
    None) the weights are computed assuming Poisson noise.  
* `instrument_response_function` :  
    The instrument response function (IRF) that is used for convolution. If no
    IRF is provided  
* `start` :  
    The start index in the IRF used for convolution. Points in the IRF before
    the start index are not used for convolution.  
* `stop` :  
    The stop index in the IRF used for convolution. Points beyond the stop index
    are not convolved.  
* `use_amplitude_threshold` :  
    If this is set to true (default value is true) the values that are smaller
    then a specified threshold are omitted  
* `amplitude_threshold` :  
    The amplitude threshold that is used if the parameter
    use_amplitude_threshold is set to true (the default value is 1e10)  
* `correct_pile_up` :  
    If this is set to true (the default value is false) the convolved model
    function is 'piled up' to match pile up artifacts in the data.  
* `excitation_period` :  
    the repetition period, .i.e, the time between subsequent excitation pulses.  
";

%feature("docstring") Decay::get_weighted_residuals "
";

%feature("docstring") Decay::evaluate "
";

%feature("docstring") Decay::get_chi2 "

Computes the chi2 for the model and the data  

The \"normal\" chi2 is the sum of the squared weighted deviations between the
data and the model.  

Parameters
----------
* `x_min` :  
    minimum index number of data / model used to compute the chi2  
* `x_max` :  
    maximum index number of data / model used to compute the chi2  
* `type` :  
    is either neyman or poisson for large count and low count data,
    respectively.  

Returns
-------
the chi2 value  
";

%feature("docstring") Decay::set "

Convenience method to update parameters that are frequently changed.  

Parameters
----------
* `irf_background` :  
* `irf_shift_channels` :  
* `areal_scatter_fraction` :  
* `constant_background` :  
* `lifetime_spectrum` :  
* `n_lifetime_spectrum` :  
";

// File: struct_l_v_double_array.xml


%feature("docstring") LVDoubleArray "
";

%feature("docstring") LVDoubleArray::~LVDoubleArray "
";

// File: struct_l_v_i32_array.xml


%feature("docstring") LVI32Array "

Structures and functions used for LabView interface.  

fit2x was originally developed as a C backend for LabView software. Therefore,
the interface with fit2x uses structures that can be accessed by Labview. In
order to make an interfacing with Python and other languages possible there is a
this files defines a set of functions that facilitate the creation of the
LabView structures.  

C++ includes: lvarrays.h
";

%feature("docstring") LVI32Array::~LVI32Array "
";

// File: struct_m_param.xml


%feature("docstring") MParam "
";

%feature("docstring") MParam::~MParam "
";

// File: classap_1_1raw__vector.xml


%feature("docstring") ap::raw_vector "
";

%feature("docstring") ap::raw_vector::raw_vector "
";

%feature("docstring") ap::raw_vector::GetData "
";

// File: classap_1_1template__1d__array.xml


%feature("docstring") ap::template_1d_array "
";

%feature("docstring") ap::template_1d_array::template_1d_array "
";

%feature("docstring") ap::template_1d_array::template_1d_array "
";

%feature("docstring") ap::template_1d_array::~template_1d_array "
";

%feature("docstring") ap::template_1d_array::setbounds "
";

%feature("docstring") ap::template_1d_array::setcontent "
";

%feature("docstring") ap::template_1d_array::getcontent "
";

%feature("docstring") ap::template_1d_array::getcontent "
";

%feature("docstring") ap::template_1d_array::getlowbound "
";

%feature("docstring") ap::template_1d_array::gethighbound "
";

%feature("docstring") ap::template_1d_array::getvector "
";

%feature("docstring") ap::template_1d_array::getvector "
";

// File: classap_1_1template__2d__array.xml


%feature("docstring") ap::template_2d_array "
";

%feature("docstring") ap::template_2d_array::template_2d_array "
";

%feature("docstring") ap::template_2d_array::template_2d_array "
";

%feature("docstring") ap::template_2d_array::~template_2d_array "
";

%feature("docstring") ap::template_2d_array::setbounds "
";

%feature("docstring") ap::template_2d_array::setcontent "
";

%feature("docstring") ap::template_2d_array::getcontent "
";

%feature("docstring") ap::template_2d_array::getcontent "
";

%feature("docstring") ap::template_2d_array::getlowbound "
";

%feature("docstring") ap::template_2d_array::gethighbound "
";

%feature("docstring") ap::template_2d_array::getcolumn "
";

%feature("docstring") ap::template_2d_array::getcolumn "
";

%feature("docstring") ap::template_2d_array::getrow "
";

%feature("docstring") ap::template_2d_array::getrow "
";

// File: namespaceap.xml

%feature("docstring") ap::abscomplex "
";

%feature("docstring") ap::conj "
";

%feature("docstring") ap::csqr "
";

%feature("docstring") ap::vdotproduct "
";

%feature("docstring") ap::vmove "
";

%feature("docstring") ap::vmove "
";

%feature("docstring") ap::vmoveneg "
";

%feature("docstring") ap::vadd "
";

%feature("docstring") ap::vadd "
";

%feature("docstring") ap::vsub "
";

%feature("docstring") ap::vsub "
";

%feature("docstring") ap::vmul "
";

%feature("docstring") ap::sign "
";

%feature("docstring") ap::randomreal "
";

%feature("docstring") ap::randominteger "
";

%feature("docstring") ap::round "
";

%feature("docstring") ap::trunc "
";

%feature("docstring") ap::ifloor "
";

%feature("docstring") ap::iceil "
";

%feature("docstring") ap::pi "
";

%feature("docstring") ap::sqr "
";

%feature("docstring") ap::maxint "
";

%feature("docstring") ap::minint "
";

%feature("docstring") ap::maxreal "
";

%feature("docstring") ap::minreal "
";

// File: namespace_phasor.xml

%feature("docstring") Phasor::compute_phasor "

Compute the phasor (g,s) for a selection of micro times  

This function computes the phasor (g,s) for a set of micro times that are
selected out of an vector. The microtimes are selected by a second vector. The
second vector speciefies which indices of the microtime vector are used to
compute the phasor.  

Parameters
----------
* `micro_times` :  
    vector of micro times  
* `idxs` :  
    vector of selected indices  
* `minimum_number_of_photons` :  
* `frequency` :  
    the frequency of the phasor  
* `g_irf` :  
    g-value of instrument response phasor  
* `s_irf` :  
    s-value of instrument response phasor  

Returns
-------
vector of length 2: first element g-value, second element s-value  
";

%feature("docstring") Phasor::compute_phasor_all "

Compute the phasor (g,s) for a all passed micro times  

Parameters
----------
* `micro_times` :  
    vector of micro times  
* `n_microtimes` :  
    number of elements in the micro time array  
* `frequency` :  
    the frequency of the phasor  

Returns
-------
vector of length 2: first element g-value, second element s-value  
";

%feature("docstring") Phasor::g "

https://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0194578.s001  

Parameters
----------
* `g_irf` :  
    g-value of instrument response phasor  
* `s_irf` :  
    s-value of instrument response phasor  
* `g_exp` :  
* `s_exp` :  

Returns
-------  
";

%feature("docstring") Phasor::s "

https://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0194578.s001  

Parameters
----------
* `g_irf` :  
* `s_irf` :  
* `g_exp` :  
* `s_exp` :  

Returns
-------  
";

// File: namespacestatistics.xml

%feature("docstring") statistics::chi2_counting "

Different chi2 measures for counting data:  

https://arxiv.org/pdf/1903.07185.pdf  

Parameters
----------
* `data` :  
* `model` :  
* `x_min` :  
* `x_max` :  
* `type` :  

Returns
-------  
";

// File: ap_8h.xml

// File: decay_8h.xml

// File: fit23_8h.xml

%feature("docstring") modelf23 "

Single exponential model function with single rotational correlation time, with
scatter contribution (BIFL scatter model)  

This function computes the fluorescence decay in the parallel and perpendicular
detection channel for a single exponential decay with a fluorescence lifetime
tau, a single rotational correlation time rho, and separate instrument response
functions for the parallel and perpendicular detection channels. The model
considers the faction of scattered light in the two detection channels by the
parameter gamma. The scattered light contribution is handled by patterns for the
light in the parallel and perpendicular detection channel.  

The instrument response function, the background, and the computed model
function are in the Jordi format, i.e., one-dimensional arrays of the parallel
and perpendicular detection channel.  

Parameters
----------
* `param` :  
    array containing the parameters of the model [0] tau, [1] gamma, [2] r0, [3]
    rho  
* `irf` :  
    instrument response function in Jordi format  
* `bg[in]` :  
    background pattern in Jordi format  
* `Nchannels[in]` :  
    number of channels (half the length of the Jordi arrays)  
* `dt[in]` :  
    time difference between two consecutive counting channels  
* `corrections[in]` :  
    [0] excitation period, [1] g factor, [2] l1, [3] l2, [4] convolution stop
    channel number  
* `mfunction[out]` :  
    output array of the computed decay in Jordi format. The output array has to
    have twice the number of channels. It needs to be allocated by beforehand.  

Returns
-------
integer (not used, 0 by default)  
";

%feature("docstring") targetf23 "

Target function (to minimize) for fit23.  

Computes the model function 23 and returns a score that quantifies the
discrepancy between the data and the model.  

Parameters
----------
* `x[in`, `out]` :  
    a vector of length that that contains the starting parameters for the
    optimization and is used to return the optimized parameters. [0]
    fluorescence lifetime - tau (in) [1] fraction of scattered light - gamma
    (in) [2] fundamental anisotropy - r0 (in) [3] rotational correlation time -
    rho (in) [4] if negative reduce contribution of background photons from
    scoring function - Soft BIFL scatter fit? (flag) (in) [5] specifies type of
    score that is returned - 2I*: P+2S? (flag), this parameter only affects the
    returned score and does not influence the fitting (in) [6] background
    corrected anisotropy - r Scatter (out) [7] anisotropy without background
    correction - r Experimental (out)  
* `pv[in]` :  
    a pointer to a MParam structure that contains the data and a set of
    corrections.  

Returns
-------
a normalized chi2  
";

%feature("docstring") fit23 "

Function optimizing parameters of model23 to the data.  

Parameters
----------
* `x[in`, `out]` :  
    a vector of length that that contains the starting parameters for the
    optimization and is used to return the optimized parameters. [0]
    fluorescence lifetime - tau (in,out) [1] fraction of scattered light - gamma
    (in,out) [2] fundamental anisotropy - r0 () [3] rotational correlation time
    - rho (in,out) [4] if negative reduce contribution of background photons
    from scoring function - Soft BIFL scatter fit? (flag) (in) [5] specifies
    type of score that is returned - 2I*: P+2S? (flag), this parameter only
    affects the returned score and does not influence the fitting (in) [6]
    background corrected anisotropy - r Scatter (out) [7] anisotropy without
    background correction - r Experimental (out)  
* `fixed` :  
    an array at least of length 4 for the parameters tau, gamma, r0, and rho
    that specifies if a parameter is optimized. If a value is set to 1, the
    parameter is optimized.  
* `p` :  
    an instance of MParam that contains all relevant information, i.e.,
    experimental data, the instrument response function, the needed corrections
    ( g-factor, l1, l2)  

Returns
-------  
";

%feature("docstring") correct_input23 "

Correct input parameters and compute anisotropy  

This function corrects the input parameters for fit23 and takes care of
unreasonable values. The fluorescence lifetime is constraint to positive values,
gamma (the fraction of scattered light) is constraint to values between 0.0 and
0.999, the rotational correlation time rho is (if the global variable fixedrho
is set to true) to the value that corresponds to the Perrin equation (for the
computed, experimental anisotropy). Moreover, this function computes the
anisotropy based on the corrected (g-factor, l1, l2, background) intensities, if
the variable return_r is set to true.  

Parameters
----------
* `x[in`, `out]` :  
    array of length 8 that contains parameters x[0] fluorescence lifetime - tau;
    x[1] fraction of scattered light - gamma; x[2] fundamental anisotropy - r0
    x[3] rotational time - rho; x[4] softbifl - flag specifying the type of bifl
    fit (not used here) x[5] p2s_twoIstar - flag specifying the type of chi2
    calculation (not used here) x[6] background corrected anisotropy x[7]
    anisotropy without background correction  
* `xm[in`, `out]` :  
    array that will contain the corrected parameters  
* `corrections[in]` :  
* `return_r[in]` :  
    if set to true (positive) computes the anisotropy and returns the scatter
    corrected and the signal (no scatter correction) anisotropy and writes the
    values to the input/output vector x.  
";

// File: fit24_8h.xml

%feature("docstring") modelf24 "

Bi-exponential model function.  

Bi-exponential model function with two fluorescence lifetimes tau1, tau2 and
amplitude of the second lifetime A2, fraction scattered light gamma, and a
constant offset. A2 (A1 + A2 = 1)  

The model function does not describe anisotropy. The decays passed as a Jordi
format are treated identical in first and the second channel of the stacked
arrays.  

mfunction[i] * (1. - gamma) / sum_m + bg[i] * gamma / sum_s + offset  

Parameters
----------
* `param` :  
    array containing the parameters of the model [0] tau1, [1] gamma, [2] tau2,
    [3] A2, [4] offset  
* `irf` :  
    instrument response function in Jordi format  
* `bg[in]` :  
    background pattern in Jordi format  
* `Nchannels[in]` :  
    number of channels (half the length of the Jordi arrays)  
* `dt[in]` :  
    time difference between two consecutive counting channels  
* `corrections[in]` :  
    [0] excitation period, [1] unused, [2] unused, [3] unused, [4] convolution
    stop channel.  
* `mfunction[out]` :  
    output array of the computed decay in Jordi format. The output array has to
    have twice the number of channels. It needs to be allocated by beforehand.  

Returns
-------  
";

%feature("docstring") targetf24 "

Target function (to minimize) for fit23.  

Parameters
----------
* `x` :  
    array containing the parameters of the model [0] tau1, [1] gamma, [2] tau2,
    [3] A2, [4] offset  
* `pv[in]` :  
    a pointer to a MParam structure that contains the data and a set of
    corrections.  

Returns
-------
a normalized chi2  
";

%feature("docstring") fit24 "

Fit a bi-exponential decay model  

This function fits a bi-exponential decay model to two decays that are stacked
using global parameters for the lifetimes and amplitudes.  

Bi-exponential model function with two fluorescence lifetimes tau1, tau2 and
amplitude of the second lifetime A2, fraction scattered light gamma, and a
constant offset. A2 (A1 + A2 = 1)  

The model function does not describe anisotropy. The decays passed as a Jordi
format are treated identical in first and the second channel of the stacked
arrays.  

The anisotropy is computed assuming that the first and the second part of the
Jordi input arrays are for parallel and perpendicular using the correction array
of the attribute p of the type MParam.  

Parameters
----------
* `x` :  
    array containing the parameters of the model [0] tau1, [1] gamma, [2] tau2,
    [3] A2, [4] offset, [5] BIFL scatter fit? (flag) - if smaller than 0 uses
    soft bifl scatter fit (seems to be unused) [6] r Scatter (output only), [7]
    r Experimental (output only)  
* `fixed` :  
    an array at least of length 5 for the parameters [0] tau1, [1] gamma, [2]
    tau2, [3] A2, [4] offset. If a value is not set to fixed the parameter is
    optimized.  
* `p` :  
    an instance of MParam that contains relevant information. Here, experimental
    data, the instrument response function, and the background decay are used.  

Returns
-------
Quality parameter 2I*  
";

%feature("docstring") correct_input24 "

Correct input parameters and compute anisotropy for fit24.  

limits (0.001 < A2 < 0.999), (0.001 < gamma < 0.999), (tau1 > 0), (tau2 > 0),
background > 0 (called offset in other places)  

Parameters
----------
* `x[in`, `out]` :  
    [0] tau1, [1] gamma [2] tau2, [3] A2, [4] background, [5] BIFL scatter fit?
    (flag, not used), [6] anisotropy r (scatter corrected, output), [7]
    anisotropy (no scatter correction, output)  
* `xm[out]` :  
    array for corrected parameters (amplied range)  
* `corrections` :  
    [1] g factor, [2] l1, [3] l3  
* `return_r` :  
    if true computes the anisotropy.  

Returns
-------  
";

// File: fit25_8h.xml

%feature("docstring") correct_input25 "

adjust parameters for fit25 and compute anisotropy  

Makes sure that (0 < gamma < 0.999) and (0<rho).  

Parameters
----------
* `x` :  
* `xm` :  
* `corrections` :  
* `return_r` :  

Returns
-------  
";

%feature("docstring") targetf25 "

Function used to compute the target value in fit 25  

This is misleadingly named target25. Fit25 selects out of a set of 4 lifetimes
the lifetime that describes best the data.  

Parameters
----------
* `x` :  
* `pv` :  

Returns
-------  
";

%feature("docstring") fit25 "

Selects the lifetime out of a set of 4 fixed lifetimes that best describes the
data.  

This function selects out of a set of 4 lifetimes tau the lifetime that fits
best the data and returns the lifetime through the parameter x[0].  

If softBIFL flag is set to (x[6] < 0) and fixed[4] is zero gamma is optimized
for each lifetime tau and the best gamma is returned by x[4]. The gamma is
fitted with fit23.  

Parameters
----------
* `x` :  
    array containing the parameters [0] tau1 output for best tau (always fixed),
    [1] tau2 (always fixed), [2] tau3 (always fixed), [3] tau4 (always fixed),
    [4] gamma (input, output), [5] fundamental anisotropy r0, [6] BIFL scatter
    fit? (flag), [7] r Scatter (output only), [8] r Experimental (output only)  
* `fixed` :  
    array that is of least of length 5. Only the element fixed[4] is used. If
    fixed[4] is zero gamma is optimized for each lifetime.  
* `p` :  
    an instance of MParam that contains all relevant information, i.e.,
    experimental data, the instrument response function, the needed corrections
    for the anisotropy (g-factor, l1, l2)  

Returns
-------  
";

// File: fit26_8h.xml

%feature("docstring") correct_input26 "

Correct input for fit 26  

Constrains the fraction x1 of the first pattern to (0 < x1 < 1).  

Parameters
----------
* `x[in]` :  
    x[0] fraction of first pattern  
* `xm[out]` :  
    xm[0] corrected fraction of first pattern  
";

%feature("docstring") targetf26 "
";

%feature("docstring") fit26 "

Pattern-fit  

Fits the fraction of a mixture of two patterns  

The two patterns are set by the attributes irf and bg of the MParam structure.  

Parameters
----------
* `x` :  
    [0] fraction of pattern 1  
* `fixed` :  
    not used  
* `p` :  
    an instance of MParam that contains the patterns. The fist pattern is
    contained in the instrument response function array, the second in the
    background, array, the experimental data is in the array expdata.  

Returns
-------  
";

// File: fits2x_8h.xml

%feature("docstring") compute_signal_and_background "

Computes the total number of photons in the parallel and perpendicular detection
channel for the background and the measured signal. The computed number of
photons are stored in the static variables Sp, Ss, Bp, Bs.  

Parameters
----------
* `p[in]` :  
    a pointer to a MParam object  
";

%feature("docstring") normM "

Normalizes the number of photons in the entire model function to the number of
experimental photons.  

Here, the Number of experimental photons is Sp + Ss (signal in parallel and
perpendicular). Sp and Ss are global variables that can be computed by
`compute_signal_and_background`.  

Parameters
----------
* `M[in`, `out]` :  
    array containing the model function in Jordi format  
* `Nchannels[in]` :  
    number of channels in the experiment(half length of M array)  
";

%feature("docstring") normM "

Normalizes a model function (that is already normalized to a unit area) to the
total number of photons in parallel and perpendicular,  

Parameters
----------
* `M[in`, `out]` :  
    array containing the model function in Jordi format  
* `s[in]` :  
    a scaling factor by which the model function is divided.  
* `Nchannels[in]` :  
    the number of channels in the model function (half length of M array)  
";

%feature("docstring") normM_p2s "

Normalizes the number of photons in the model function for Ss and Sp
individually to the number of experimental photons in Ss and Sp.  

Here, the number of experimental photons are global variables that can be
computed by `compute_signal_and_background`.  

Parameters
----------
* `M` :  
    array[in,out] containing the model function in Jordi format  
* `Nchannels[in]` :  
    number of channels in the experiment (half length of M array)  
";

// File: fsconv_8h.xml

%feature("docstring") rescale "

Convolution, scaling, and lamp shift routines.  

Scale model function to the data (old version)  

This function rescales the model function (fit) to the data by the number of
photons between a start and a stop micro time counting channel. The number of
photons between start and stop are counted and the model function is scaled to
match the data by area.  

This rescaling function does not consider the noise in the data when rescaling
the model.  

Parameters
----------
* `fit[in`, `out]` :  
    model function that is scaled (modified in-place)  
* `decay[in]` :  
    the experimental data to which the model function is scaled  
* `scale[out]` :  
    the scaling parameter (the factor) by which the model function is
    multiplied.  
* `start[in]` :  
    The start micro time channel  
* `stop[in]` :  
    The stop micro time channel  
";

%feature("docstring") rescale_w "

Scale model function to the data (with weights)  

This function rescales the model function (fit) to the data by the number of
photons between a start and a stop micro time counting channel. The number of
photons between start and stop are counted and the model function is scaled to
match the data by area considering the noise of the data.  

The scaling factor is computed by:  

scale = sum(fit*decay/w^2)/sum(fit^2/w^2)  

Parameters
----------
* `fit[in`, `out]` :  
    model function that is scaled (modified in-place)  
* `decay[in]` :  
    the experimental data to which the model function is scaled  
* `w_sq[in]` :  
    squared weights of the data.  
* `scale[out]` :  
    the scaling parameter (the factor) by which the model function is
    multiplied.  
* `start[in]` :  
    The start micro time channel  
* `stop[in]` :  
    The stop micro time channel  
";

%feature("docstring") rescale_w_bg "

Scale model function to the data (with weights and background)  

This function scales the model function (fit) to the data by the number of
photons between a start and a stop micro time counting channel. The number of
photons between start and stop are counted and the model function is scaled to
match the data by area considering the noise of the data and a constant offset
of the data.  

scale = sum(fit*(decay-bg)/w^2)/sum(fit^2/w^2)  

Parameters
----------
* `fit[in`, `out]` :  
    model function that is scaled (modified in-place)  
* `decay[in]` :  
    the experimental data to which the model function is scaled  
* `w_sq[in]` :  
    squared weights of the data.  
* `bg[in]` :  
    constant background of the data  
* `scale[out]` :  
    the scaling parameter (the factor) by which the model function is
    multiplied.  
* `start[in]` :  
    The start micro time channel  
* `stop[in]` :  
    The stop micro time channel  
";

%feature("docstring") fconv "

Convolve lifetime spectrum with instrument response (fast convolution, low
repetition rate)  

This function computes the convolution of a lifetime spectrum (a set of
lifetimes with corresponding amplitudes) with a instrument response function
(irf). This function does not consider periodic excitation and is suited for
experiments at low repetition rate.  

Parameters
----------
* `fit[out]` :  
    model function. The convoluted decay is written to this array  
* `x[in]` :  
    lifetime spectrum (amplitude1, lifetime1, amplitude2, lifetime2, ...)  
* `lamp[in]` :  
    instrument response function  
* `numexp[in]` :  
    number of fluorescence lifetimes  
* `start[in]` :  
    start micro time index for convolution (not used)  
* `stop[in]` :  
    stop micro time index for convolution.  
* `dt[in]` :  
    time difference between two micro time channels  
";

%feature("docstring") fconv_per "

Convolve lifetime spectrum with instrument response (fast convolution, high
repetition rate)  

This function computes the convolution of a lifetime spectrum (a set of
lifetimes with corresponding amplitudes) with a instrument response function
(irf). This function does consider periodic excitation and is suited for
experiments at high repetition rate.  

Parameters
----------
* `fit[out]` :  
    model function. The convoluted decay is written to this array  
* `x[in]` :  
    lifetime spectrum (amplitude1, lifetime1, amplitude2, lifetime2, ...)  
* `lamp[in]` :  
    instrument response function  
* `numexp[in]` :  
    number of fluorescence lifetimes  
* `start[in]` :  
    start micro time index for convolution (not used)  
* `stop[in]` :  
    stop micro time index for convolution.  
* `n_points` :  
    number of points in the model function.  
* `period` :  
    excitation period in units of the fluorescence lifetimes (typically
    nanoseconds)  
* `dt[in]` :  
    time difference between two micro time channels  
";

%feature("docstring") fconv_per_cs "

Convolve lifetime spectrum - fast convolution, high repetition rate, with
convolution stop.  

fast convolution, high repetition rate, with convolution stop for Paris  

Parameters
----------
* `fit[out]` :  
    model function. The convoluted decay is written to this array  
* `x[in]` :  
    lifetime spectrum (amplitude1, lifetime1, amplitude2, lifetime2, ...)  
* `lamp[in]` :  
    instrument response function  
* `numexp[in]` :  
    number of fluorescence lifetimes  
* `stop[in]` :  
    stop micro time index for convolution.  
* `n_points` :  
    number of points in the model function.  
* `period` :  
    excitation period in units of the fluorescence lifetimes (typically
    nanoseconds)  
* `conv_stop` :  
    convolution stop micro channel number  
* `dt[in]` :  
    time difference between two micro time channels  
";

%feature("docstring") fconv_ref "

Convolve lifetime spectrum - fast convolution with reference compound decay.  

This function convolves a set of fluorescence lifetimes and with associated
amplitudes with an instrument response function. The provided amplitudes are
scaled prior to the convolution by area using a reference fluorescence lifetime.
The amplitudes are computed by  

amplitude_corrected = a * ( 1 /tauref - 1 / tau)  

where a and tau are provided amplitudes.  

Parameters
----------
* `fit[out]` :  
    model function. The convoluted decay is written to this array  
* `x[in]` :  
    lifetime spectrum (amplitude1, lifetime1, amplitude2, lifetime2, ...)  
* `lamp[in]` :  
    instrument response function  
* `numexp[in]` :  
    number of fluorescence lifetimes  
* `start[in]` :  
    start micro time index for convolution (not used)  
* `stop[in]` :  
    stop micro time index for convolution.  
* `tauref` :  
    a reference lifetime used to rescale the amplitudes of the fluorescence
    lifetime spectrum  
* `dt[in]` :  
    time difference between two micro time channels  
";

%feature("docstring") sconv "

Convolve fluorescence decay curve with irf - slow convolution.  

This function computes a convolved model function for a fluorescence decay
curve.  

Parameters
----------
* `fit` :  
    convolved model function  
* `p` :  
    model function before convolution - fluorescence decay curve  
* `lamp` :  
    instrument response function  
* `start` :  
    start index of the convolution  
* `stop` :  
    stop index of the convolution  
";

%feature("docstring") shift_lamp "

shift instrumnet response function  

Parameters
----------
* `lampsh` :  
* `lamp` :  
* `ts` :  
* `n_points` :  
* `out_value` :  
    the value of the shifted response function outside of the valid indices  
";

%feature("docstring") add_pile_up "

Correct the model function for pile up.  

Add pile up to a model function. The pile-up model follows the description by
Coates, 1968, eq. 2  

p = data / (n_excitation_pulses - np.cumsum(data)) Coates, 1968, eq. 4  

Reference: Coates, P.: The correction for photonpile-up’ in the measurement of
radiative lifetimes. J. Phys. E: Sci. Instrum. 1(8), 878–879 (1968)  

Parameters
----------
* `model[in`, `out]` :  
    The array containing the model function  
* `n_model[in]` :  
    Number of elements in the model array  
* `data[in]` :  
    The array containing the experimental decay  
* `n_data[in]` :  
    number of elements in experimental decay  
* `repetition_rate[in]` :  
    The repetition-rate in MHz  
* `dead_time[in]` :  
    The dead-time of the detection system in nanoseconds  
* `measurement_time[in]` :  
    The measurement time in seconds  
* `pile_up_model[in]` :  
    The model used to compute the pile up distortion of the data (currently only
    Coates)  
";

%feature("docstring") fconv_per_cs_time_axis "

Compute the fluorescence decay for a lifetime spectrum and a instrument response
function considering periodic excitation.  

Fills the pre-allocated output array `output_decay` with a fluorescence
intensity decay defined by a set of fluorescence lifetimes defined by the
parameter `lifetime_spectrum`. The fluorescence decay will be convolved (non-
periodically) with an instrumental response function that is defined by
`instrument_response_function`.  

This function calculates a fluorescence intensity model_decay that is convolved
with an instrument response function (IRF). The fluorescence intensity
model_decay is specified by its fluorescence lifetime spectrum, i.e., an
interleaved array containing fluorescence lifetimes with corresponding
amplitudes.  

This convolution only works with evenly linear spaced time axes.  

Parameters
----------
* `inplace_output[in`, `out]` :  
    Inplace output array that is filled with the values of the computed
    fluorescence intensity decay model  
* `n_output[in]` :  
    Number of elements in the output array  
* `time_axis[in]` :  
    the time-axis of the model_decay  
* `n_time_axis[in]` :  
    length of the time axis  
* `irf[in]` :  
    the instrument response function array  
* `n_irf[in]` :  
    length of the instrument response function array  
* `lifetime_spectrum[in]` :  
    Interleaved array of amplitudes and fluorescence lifetimes of the form
    (amplitude, lifetime, amplitude, lifetime, ...)  
* `n_lifetime_spectrum[in]` :  
    number of elements in the lifetime spectrum  
* `convolution_start[in]` :  
    Start channel of convolution (position in array of IRF)  
* `convolution_stop[in]` :  
    convolution stop channel (the index on the time-axis)  
* `use_amplitude_threshold[in]` :  
    If this value is True (default False) fluorescence lifetimes in the lifetime
    spectrum which have an amplitude with an absolute value of that is smaller
    than `amplitude_threshold` are not omitted in the convolution.  
* `amplitude_threshold[in]` :  
    Threshold value for the amplitudes  
* `period` :  
    Period of repetition in units of the lifetime (usually, nano-seconds)  
";

%feature("docstring") fconv_cs_time_axis "

Compute the fluorescence decay for a lifetime spectrum and a instrument response
function.  

Fills the pre-allocated output array `output_decay` with a fluorescence
intensity decay defined by a set of fluorescence lifetimes defined by the
parameter `lifetime_spectrum`. The fluorescence decay will be convolved (non-
periodically) with an instrumental response function that is defined by
`instrument_response_function`.  

This function calculates a fluorescence intensity model_decay that is convolved
with an instrument response function (IRF). The fluorescence intensity
model_decay is specified by its fluorescence lifetime spectrum, i.e., an
interleaved array containing fluorescence lifetimes with corresponding
amplitudes.  

This convolution works also with uneven spaced time axes.  

Parameters
----------
* `inplace_output[in`, `out]` :  
    Inplace output array that is filled with the values of the computed
    fluorescence intensity decay model  
* `n_output[in]` :  
    Number of elements in the output array  
* `time_axis[in]` :  
    the time-axis of the model_decay  
* `n_time_axis[in]` :  
    length of the time axis  
* `irf[in]` :  
    the instrument response function array  
* `n_irf[in]` :  
    length of the instrument response function array  
* `lifetime_spectrum[in]` :  
    Interleaved array of amplitudes and fluorescence lifetimes of the form
    (amplitude, lifetime, amplitude, lifetime, ...)  
* `n_lifetime_spectrum[in]` :  
    number of elements in the lifetime spectrum  
* `convolution_start[in]` :  
    Start channel of convolution (position in array of IRF)  
* `convolution_stop[in]` :  
    convolution stop channel (the index on the time-axis)  
* `use_amplitude_threshold[in]` :  
    If this value is True (default False) fluorescence lifetimes in the lifetime
    spectrum which have an amplitude with an absolute value of that is smaller
    than `amplitude_threshold` are not omitted in the convolution.  
* `amplitude_threshold[in]` :  
    Threshold value for the amplitudes  
";

// File: i__lbfgs_8h.xml

%feature("docstring") fjac1 "
";

%feature("docstring") fgrad1 "
";

%feature("docstring") fjac2 "
";

%feature("docstring") fgrad2 "
";

%feature("docstring") fjac4 "
";

%feature("docstring") fgrad4 "
";

// File: lbfgs_8h.xml

%feature("docstring") funcgrad "
";

%feature("docstring") lbfgsminimize "
";

// File: lvarrays_8h.xml

%feature("docstring") CreateLVI32Array "

Parameters
----------
* `len` :  

Returns
-------  
";

%feature("docstring") CreateLVDoubleArray "

Parameters
----------
* `len` :  

Returns
-------  
";

%feature("docstring") CreateMParam "
";

// File: phasor_8h.xml

// File: statistics_8h.xml

// File: two_istar_8h.xml

%feature("docstring") init_fact "

Initialize an array containing pre-computed logratithms  
";

%feature("docstring") loggammaf "

Approximation of log(gamma function). See wikipedia  

https://en.wikipedia.org/wiki/Gamma_function#The_log-gamma_function  

Parameters
----------
* `t` :  
    input of the gamma function  

Returns
-------
approximation of the logarithm of the gamma function  
";

%feature("docstring") wcm "

log-likelihood w(C|m) for Cp + 2Cs  

Parameters
----------
* `C` :  
    number of counts in channel  
* `m` :  
    model function  

Returns
-------
log-likelihood w(C|m) for Cp + 2Cs  
";

%feature("docstring") wcm_p2s "

Compute the -log-likelihood for Cp + 2Cs of a single micro time channel.  

Compute score of model counts in a parallel and perpendicular detection channel
and the experimental counts for a micro time channel.  

This function computes a score for the experimental counts (C) in a channel
where the experimental counts were computed by the sum of the counts in the
parallel (P) and the perpendicular (S) channel by the equation C = P + 2 S.  

This function considers that the number of counts C = P + 2S is not Poissonian.
The score relates to a maximum likelihood function.  

Parameters
----------
* `C` :  
    number of experimental counts (P + 2 S) in a micro time channel  
* `mp` :  
    number of counts of the model in parallel detection channel  
* `ms` :  
    number of counts of the model in the perpendicular detection channel  

Returns
-------  
";

%feature("docstring") Wcm_p2s "

Compute the overall -log-likelihood for Cp + 2Cs for all micro time channels  

Parameters
----------
* `C` :  
    array of experimental counts in Jordi format  
* `M` :  
    array model function in Jordi format  
* `Nchannels` :  
    number of micro time channels in parallel and perpendicular (half the number
    of elements in C and M).  

Returns
-------
-log-likelihood for Cp + 2Cs for all micro time channels  
";

%feature("docstring") twoIstar_p2s "

Compute overall 2I* for Cp + 2Cs  

This function computes the overall 2I* for the model function Cp + 2Cs that is
computed by parallel signal (Cp) and the perpendicular signal (Cs). For the
definition of function 2I* see \"An Experimental Comparison of the Maximum
Likelihood Estimation and Nonlinear Least-Squares Fluorescence Lifetime Analysis
of Single Molecules, Michael Maus, Mircea Cotlet, Johan Hofkens, Thomas Gensch,
Frans C. De Schryver, J. Schaffer, and C. A. M. Seidel, Anal. Chem. 2001, 73, 9,
2078–2086\".  

Parameters
----------
* `C` :  
    array of experimental counts in Jordi format  
* `M` :  
    array model function in Jordi format  
* `Nchannels` :  
    number of micro time channels in parallel and perpendicular (half the number
    of elements in C and M).  

Returns
-------
2I* for Cp + 2Cs  
";

%feature("docstring") twoIstar "

Compute overall 2I* for Cp & Cs  

This function computes 2I* for Cp and Cs. Cp and Cs are the model signals in the
parallel and perpendicular channel. Contrary to twoIstar_p2s the overall 2I* is
the sum of 2I* for Cp and Cs.  

Parameters
----------
* `C` :  
    array of experimental counts in Jordi format  
* `M` :  
    array model function in Jordi format  
* `Nchannels` :  
    number of micro time channels in parallel and perpendicular (half the number
    of elements in C and M).  

Returns
-------
2I* for Cp & Cs  
";

%feature("docstring") Wcm "

Compute overall -log-likelihood for Cp & Cs  

Parameters
----------
* `C` :  
    array of experimental counts in Jordi format  
* `M` :  
    array model function in Jordi format  
* `Nchannels` :  
    number of micro time channels in parallel and perpendicular (half the number
    of elements in C and M).  

Returns
-------
-log-likelihood for Cp & Cs  
";

// File: dir_d44c64559bbebec7f509842c48db8b23.xml


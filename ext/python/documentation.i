
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

// File: struct_l_v_double_array.xml


%feature("docstring") LVDoubleArray "
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

// File: struct_m_param.xml


%feature("docstring") MParam "
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

// File: ap_8h.xml

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
    a pointer to a MParam structure.  

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

Bi-exponential model function  

Parameters
----------
* `param` :  
* `irf` :  
* `bg` :  
* `Nchannels` :  
* `dt` :  
* `corrections` :  
* `mfunction` :  

Returns
-------  
";

%feature("docstring") targetf24 "

Parameters
----------
* `x` :  
* `pv` :  

Returns
-------  
";

%feature("docstring") fit24 "

Parameters
----------
* `x` :  
* `fixed` :  
* `p` :  

Returns
-------  
";

%feature("docstring") correct_input24 "

Parameters
----------
* `x` :  
* `xm` :  
* `corrections` :  
* `return_r` :  

Returns
-------  
";

// File: fit25_8h.xml

%feature("docstring") correct_input25 "

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

Parameters
----------
* `x` :  
* `pv` :  

Returns
-------  
";

%feature("docstring") fit25 "

Parameters
----------
* `x` :  
* `fixed` :  
* `p` :  

Returns
-------  
";

// File: fit26_8h.xml

%feature("docstring") correct_input26 "
";

%feature("docstring") targetf26 "
";

%feature("docstring") fit26 "
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

rescaling -- old version. sum(fit)->sum(decay)  

Parameters
----------
* `fit` :  
* `decay` :  
* `scale` :  
* `start` :  
* `stop` :  
";

%feature("docstring") rescale_w "

rescaling new  

Parameters
----------
* `fit` :  
* `decay` :  
* `w_sq` :  
* `scale` :  
* `start` :  
* `stop` :  
";

%feature("docstring") rescale_w_bg "

rescaling -- new version + background. scale = sum(fit*decay/w^2)/sum(fit^2/w^2)  

Parameters
----------
* `fit` :  
* `decay` :  
* `w_sq` :  
* `bg` :  
* `scale` :  
* `start` :  
* `stop` :  
";

%feature("docstring") fconv "

fast convolution  

Parameters
----------
* `fit` :  
* `x` :  
* `lamp` :  
* `numexp` :  
* `start` :  
* `stop` :  
";

%feature("docstring") fconv_per "

fast convolution, high repetition rate  

Parameters
----------
* `fit` :  
* `x` :  
* `lamp` :  
* `numexp` :  
* `start` :  
* `stop` :  
* `n_points` :  
* `period` :  
";

%feature("docstring") fconv_per_cs "

fast convolution, high repetition rate, with convolution stop  

Parameters
----------
* `fit` :  
* `x` :  
* `lamp` :  
* `numexp` :  
* `stop` :  
* `n_points` :  
* `period` :  
* `conv_stop` :  
";

%feature("docstring") fconv_ref "

fast convolution with reference compound decay  

Parameters
----------
* `fit` :  
* `x` :  
* `lamp` :  
* `numexp` :  
* `start` :  
* `stop` :  
* `tauref` :  
";

%feature("docstring") sconv "

slow convolution  

Parameters
----------
* `fit` :  
* `p` :  
* `lamp` :  
* `start` :  
* `stop` :  
";

%feature("docstring") shift_lamp "

shifting lamp  

Parameters
----------
* `lampsh` :  
* `lamp` :  
* `ts` :  
* `n_points` :  
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


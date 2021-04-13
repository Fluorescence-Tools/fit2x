#include "DecayConvolution.h"


void DecayConvolution::convolve_lifetimes(){
#ifdef VERBOSE_FIT2X
    std::clog << "DecayConvolution::convolve_lifetimes" << std::endl;
#endif
    DecayCurve irfc = *get_corrected_irf();
    auto ls = get_lifetime_spectrum();
    double* lt; int ln;
    ls->get_lifetime_spectrum(&lt, &ln);
    int cv_start = get_convolution_start();
    int cv_stop = get_convolution_stop();
    double dt = decay.get_dx();
    int nl = ln / 2;

    decay.resize(irfc.size());
#ifdef VERBOSE_FIT2X
    std::clog << "-- convolution_method:" << convolution_method << std::endl;
    std::clog << "-- convolution start: " << cv_start << std::endl;
    std::clog << "-- convolution stop: " << cv_stop << std::endl;
    std::clog << "-- dt: " << dt << std::endl;
    std::clog << "-- excitation_period: " << excitation_period << std::endl;
    std::clog << "-- number of lifetimes: " << ln / 2<< std::endl;
    std::clog << "-- lifetime spectrum: ";
    for(int i=0; i< ln; i++) std::clog << lt[i] << ", ";
    std::clog << std::endl;
    std::clog << "-- decay.size(): " << decay.size() << std::endl;
    std::clog << "-- irfc.size(): " << irfc.size() << std::endl;
    std::clog << "-- Corrected irf: ";
    for(auto &v: irfc.y) std::clog  << v << " ";
    std::clog  << std::endl;
#endif

    if(convolution_method == ConvFastPeriodicTime){
        fconv_per_cs_time_axis(
                decay.y.data(), decay.y.size(),
                irfc.x.data(), irfc.x.size(),
                irfc.y.data(), irfc.y.size(),
                lt, ln,
                cv_start, cv_stop,
                excitation_period
        );
    } else if(convolution_method == ConvFastTime){
        fconv_cs_time_axis(
                decay.y.data(), decay.y.size(),
                decay.x.data(), decay.x.size(),
                irfc.y.data(), irfc.y.size(),
                lt, ln,
                cv_start, cv_stop
        );
    } else {
        if (convolution_method == ConvFastPeriodic) {
            fconv_per(
                    decay.y.data(), lt,
                    irfc.y.data(), nl ,
                    cv_start, cv_stop,
                    decay.y.size(), excitation_period,
                    dt
            );
        } else if (convolution_method == ConvFast) {
            std::cout << "fconv" << std::endl;
            fconv(
                    decay.y.data(), lt,
                    irfc.y.data(), nl,
                    cv_start, cv_stop,
                    dt
            );
        } else if (convolution_method == ConvFastAVX) {
            fconv_avx(
                    decay.y.data(), lt,
                    irfc.y.data(), nl,
                    cv_start, cv_stop,
                    dt
            );
        }  else if (convolution_method == ConvFastPerAVX) {
            fconv_per_avx(
                    decay.y.data(), lt,
                    irfc.y.data(), nl,
                    cv_start, cv_stop,
                    decay.y.size(), excitation_period,
                    dt
            );
        }
    }
#ifdef VERBOSE_FIT2X
    std::clog << "-- decay: ";
    for(auto &v: decay.y) std::clog << " " << v;
    std::clog << std::endl;
#endif
}


void DecayConvolution::add_scatter(){
    DecayCurve sc = use_corrected_irf_as_scatter ?
                    *get_corrected_irf() : *get_irf() ;
    int n; double* v;
    DecayCurve::add_arrays(
            &v, &n,
            decay.y.data(), decay.y.size(),
            sc.y.data(), sc.y.size(),
            scatter_fraction,
            convolution_start, convolution_stop
    );
}

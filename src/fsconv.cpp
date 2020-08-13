// fast and slow convolution routines, autoscaling and lamp shift
// Version 2006.03.26

#include "fsconv.h"
#include "fits2x.h"

/* rescaling -- old version. sum(fit)->sum(decay) */
void rescale(double *fit, double *decay, double *scale, int start, int stop)
{
    double sumfit=0., sumcurve=0.;
    /* scaling */
    if (*scale==0.){
      for (int i=start; i<=stop; i++){
        sumfit += fit[i];
        sumcurve += decay[i];
      }
      if (sumfit!=0.) *scale=sumcurve/sumfit;
     }
#if VERBOSE
    std::cout << "RESCALE" << std::endl;
    std::cout << "start / stop: " << start << " / " << stop << std::endl;
    std::cout << "sumfit: " << sumfit << std::endl;
    std::cout << "sumcurve: " << sumcurve << std::endl;
    std::cout << "scale: " << *scale << std::endl;
#endif
    for (int i=start; i<=stop; i++)
      fit[i] *= *scale;
}

/* rescaling -- new version. scale = sum(fit*decay/w^2)/sum(fit^2/w^2) */
void rescale_w(double *fit, double *decay, double *w_sq, double *scale, int start, int stop)
{
    int i;
    double sumnom=0., sumdenom=0.;

    /* scaling */
    if (*scale==0.){
      for (i=start; i<=stop; i++){
        if (decay[i]!=0.) {
          sumnom += fit[i]*decay[i]/w_sq[i];
          sumdenom += fit[i]*fit[i]/w_sq[i];
        }
      }
      if (sumdenom!=0.) *scale = sumnom/sumdenom;
     }
    for (i=start; i<=stop; i++)
      fit[i] *= *scale;

}

/* rescaling -- new version + background. scale = sum(fit*decay/w^2)/sum(fit^2/w^2) */
void rescale_w_bg(double *fit, double *decay, double *w_sq, double bg, double *scale, int start, int stop)
{
#if VERBOSE
    std::clog << "RESCALE_W_BG" << std::endl;
    std::clog << "-- initial scale: " << *scale << std::endl;
    std::clog << "w_sq [:64]: ";
    for(int i=0; i<64; i++) std::clog << w_sq[i] << " ";
    std::clog << std::endl;
    std::clog << "decay [:64]: ";
    for(int i=0; i<64; i++) std::clog << decay[i] << " ";
    std::clog << std::endl;
    std::clog << "fit [:64]: ";
    for(int i=0; i<64; i++) std::clog << fit[i] << " ";
    std::clog << std::endl;
#endif
    int i;
    double sumnom=0., sumdenom=0.;
    /* scaling */
    if (*scale==0.){
      for (i=start; i<=stop; i++){
        if (decay[i]!=0.) {
          sumnom += fit[i]*(decay[i]-bg)/w_sq[i];
          sumdenom += fit[i]*fit[i]/w_sq[i];
        }
      }
      if (sumdenom!=0.) *scale = sumnom/sumdenom;
     }
#if VERBOSE
    std::clog << "-- final scale: " << *scale << std::endl;
#endif
    for (i=start; i<=stop; i++)
      fit[i] *= *scale;
}

/* fast convolution */
void fconv(double *fit, double *x, double *lamp, int numexp, int start, int stop, double dt) {
#if VERBOSE
    std::cout << "FCONV" << std::endl;
#endif
    int ne, i;
    double fitcurr, expcurr, deltathalf = dt * 0.5;
    for (i = 0; i <= stop; i++) fit[i] = 0;
    /* convolution */
    for (ne = 0; ne < numexp; ne++) {
#if VERBOSE
        std::cout << "-- tau: " << x[2 * ne + 1] << std::endl;
        std::cout << "-- x: " << x[2 * ne] << std::endl;
#endif
        expcurr = exp(-dt / x[2 * ne + 1]);
        fitcurr = 0;
        for (i = 1; i <= stop; i++) {
            fitcurr = (fitcurr + deltathalf * lamp[i - 1]) * expcurr + deltathalf * lamp[i];
            fit[i] += fitcurr * x[2 * ne];
        }
    }

}

/* fast convolution, high repetition rate */
void fconv_per(double *fit, double *x, double *lamp, int numexp, int start, int stop,
           int n_points, double period, double dt)
{
    int ne, i, lamp_start = 0, 
      stop1, period_n = (int)ceil(period/dt-0.5);
    double fitcurr, expcurr, tail_a, deltathalf = dt*0.5;

    while (lamp[lamp_start++]==0);
    for (i=0;i<=stop;i++) fit[i]=0;

    stop1 = (period_n+lamp_start > n_points-1) ? n_points-1 : period_n+lamp_start;

    /* convolution */
    for (ne=0; ne<numexp; ne++) {
      expcurr = exp(-dt/x[2*ne+1]);
      tail_a = 1./(1.-exp(-period/x[2*ne+1]));
      fitcurr = 0;
      for (i=1; i<=stop1; i++){
        fitcurr=(fitcurr + deltathalf*lamp[i-1])*expcurr + deltathalf*lamp[i];
        fit[i] += fitcurr*x[2*ne];
      }
      fitcurr *= exp(-(period_n - stop1 + start)*dt/x[2*ne+1]);
      for (i=start; i<=stop; i++){
        fitcurr *= expcurr;
        fit[i] += fitcurr*x[2*ne]*tail_a;
      }
    }

}

/* fast convolution, high repetition rate, with convolution stop for Paris */
void fconv_per_cs(double *fit, double *x, double *lamp, int numexp, int stop,
           int n_points, double period, int conv_stop, double dt)
{
    int ne, i, 
      stop1, period_n = (int)ceil(period/dt-0.5);
    double fitcurr, expcurr, tail_a, deltathalf = dt*0.5;

    for (i=0; i<=stop; i++) fit[i]=0;
    stop1 = (period_n > n_points-1) ? n_points-1 : period_n;

    /* convolution */
    for (ne=0; ne<numexp; ne++) {
      expcurr = exp(-dt/x[2*ne+1]);
      tail_a = 1./(1.-exp(-period/x[2*ne+1]));
      fitcurr = 0.;
      fit[0] += deltathalf*lamp[0]*(expcurr + 1.)*x[2*ne];
      for (i=1; i<=conv_stop; i++) {
        fitcurr=(fitcurr + deltathalf*lamp[i-1])*expcurr + deltathalf*lamp[i];
        fit[i] += fitcurr*x[2*ne];
      }
      for (; i<=stop1; i++) {
        fitcurr=fitcurr*expcurr;
        fit[i] += fitcurr*x[2*ne];
      }
      fitcurr *= exp(-(period_n - stop1)*dt/x[2*ne+1]);
      for (i=0; i<=stop; i++) {
        fitcurr *= expcurr;
        fit[i] += fitcurr*x[2*ne]*tail_a;
      }
    }
}

/* fast convolution with reference compound decay */
void fconv_ref(double *fit, double *x, double *lamp, int numexp, int start, int stop, double tauref, double dt)
{
    int ne,i;
    double fitcurr, expcurr, deltathalf = dt*0.5, correct_a, sum_a=0;

    for (i=0;i<=stop;i++) fit[i]=0; 

    /* convolution */
    for (ne=0; ne<numexp; ne++) {
      expcurr = exp(-dt/x[2*ne+1]);
      correct_a = x[2*ne]*(1/tauref-1/x[2*ne+1]);
      sum_a += x[2*ne];
      fitcurr = 0;
      for (i=1; i<=stop; i++){
        fitcurr=(fitcurr + deltathalf*lamp[i-1])*expcurr + deltathalf*lamp[i];
        fit[i] += fitcurr * correct_a;
      }
    }
    for (i=1; i<=stop; i++) fit[i] += lamp[i] * sum_a;

}

/* slow convolution */
void sconv(double *fit, double *p, double *lamp, int start, int stop)
{
    int i,j;

    /* convolution */
    for (i=start; i<=stop; i++){
       fit[i] = 0.5 * lamp[0] * p[i];
       for (j=1; j<i; j++) fit[i] += lamp[j] * p[i-j];
       fit[i] += 0.5 * lamp[i] * p[0];
       fit[i]=fit[i];
    }
    fit[0] = 0;

}


/* shifting lamp */
void shift_lamp(double *lampsh, double *lamp, double ts, int n_points, double out_value)
{
    int tsint = (int)(floor(ts));
    double tsdbl = ts-(double)tsint;
    int out_left=0, out_right=0, j; 

    if (tsint<0) out_left = -tsint;
    if (tsint+1>0) out_right = tsint+1;

    for(j=0; j<out_left; j++) lampsh[j]=out_value;
    for(j=out_left; j<(n_points-out_right); j++) 
       lampsh[j]=lamp[j+tsint]*(1-tsdbl)+lamp[j+tsint+1]*(tsdbl);
    for(j=(n_points-out_right); j<n_points; j++) lampsh[j]=out_value;

}


void add_pile_up_to_model(
        double* model, int n_model,
        double* data, int n_data,
        double repetition_rate,
        double dead_time,
        double measurement_time,
        std::string pile_up_model
){
#if VERBOSE
    std::clog << "ADD PILE-UP" << std::endl;
    std::clog << "-- Repetition_rate [MHz]: " << repetition_rate << std::endl;
    std::clog << "-- Dead_time [ns]: " << dead_time << std::endl;
    std::clog << "-- Measurement_time [s]: " << measurement_time << std::endl;
    std::clog << "-- n_data: " << n_data << std::endl;
    std::clog << "-- n_model: " << n_model << std::endl;
#endif
    if(pile_up_model == "coates"){
#if VERBOSE
        std::clog << "-- pile_up_model: " << pile_up_model << std::endl;
#endif
        repetition_rate *= 1e6;
        dead_time *= 1e-9;
        std::vector<double> cum_sum(n_data);
        std::partial_sum(data, data + n_data, cum_sum.begin(), std::plus<double>());
        long n_pulse_detected = cum_sum[cum_sum.size() - 1];
        double total_dead_time = n_pulse_detected * dead_time;
        double live_time = measurement_time - total_dead_time;
        long n_excitation_pulses = std::max(live_time * repetition_rate, (double) n_pulse_detected);
#if VERBOSE
        std::clog << "-- live_time [s]: " << live_time << std::endl;
        std::clog << "-- total_dead_time [s]: " << total_dead_time << std::endl;
        std::clog << "-- n_pulse_detected [#]: " << n_pulse_detected << std::endl;
        std::clog << "-- n_excitation_pulses [#]: " << n_excitation_pulses << std::endl;
#endif
        // Coates, 1968, eq. 2 & 4
        std::vector<double> rescaled_data(n_data);
        for(int i=0;i<n_data;i++)
            rescaled_data[i] = -std::log(1.0 - data[i] / (n_excitation_pulses - cum_sum[i]));
        for(int i=0;i<n_data;i++)
            rescaled_data[i] = (rescaled_data[i] == 0) ? 1.0 : rescaled_data[i];
        // rescale model function to preserve data counting statistics
        std::vector<double> sf(n_data);
        for(int i=0;i<n_data;i++)
            sf[i] = data[i] / rescaled_data[i];
        double s = std::accumulate(sf.begin(),sf.end(),0.0);
        for(int i=0;i<n_data;i++)
            model[i] = model[i] * (sf[i] / s * n_data);
    }
}


void discriminate_small_amplitudes(
        double* lifetime_spectrum, int n_lifetime_spectrum,
        double amplitude_threshold
        ){
    int number_of_exponentials = n_lifetime_spectrum / 2;
#if VERBOSE
    std::clog << "APPLY_AMPLITUDE_THRESHOLD" << std::endl;
    std::clog << "-- amplitude_threshold spectrum: " << amplitude_threshold << std::endl;
    std::clog << "-- lifetime spectrum before: ";
    for (int i=0; i < number_of_exponentials * 2; i++){
        std::clog << lifetime_spectrum[i] << ' ';
    }
    std::clog << std::endl;
#endif
    for(int ne = 0; ne<number_of_exponentials; ne++){
        double amplitude = lifetime_spectrum[2 * ne];
        if(std::abs(amplitude) < amplitude_threshold){
            lifetime_spectrum[2 * ne] = 0.0;
        }
    }
#if VERBOSE
    std::clog << "-- lifetime spectrum after: ";
    for (int i=0; i < number_of_exponentials * 2; i++){
        std::clog << lifetime_spectrum[i] << ' ';
    }
    std::clog << std::endl;
#endif
}


/* fast convolution, high repetition rate, with time axis */
void fconv_per_cs_time_axis(
        double* model, int n_model,
        double* time_axis, int n_time_axis,
        double *instrument_response_function, int n_instrument_response_function,
        double* lifetime_spectrum, int n_lifetime_spectrum,
        int convolution_start,
        int convolution_stop,
        double period
){
    convolution_stop = convolution_stop > 0 ?
            std::min({n_time_axis, n_instrument_response_function, n_model, convolution_stop}):
            std::min({n_time_axis, n_instrument_response_function, n_model});
    int number_of_exponentials = n_lifetime_spectrum / 2;
    double dt = time_axis[1] - time_axis[0];
    double dt_2 = dt / 2;
    int period_n = std::ceil(period / dt - 0.5);
    int stop1 = std::min(convolution_stop, period_n);
    for(int i=0; i < n_model; i++) model[i] = 0;
    convolution_start = std::max(convolution_start, 1);
#if VERBOSE
    std::clog << "convolve_lifetime_spectrum_periodic..." << std::endl;
    std::clog << "-- number_of_exponentials: " << number_of_exponentials << std::endl;
    std::clog << "-- convolution_start: " << convolution_start << std::endl;
    std::clog << "-- convolution_stop: " << convolution_stop << std::endl;
    std::clog << "-- period: " << period << std::endl;
    std::clog << "-- n_model: " << n_model << std::endl;
#endif
    for(int ne=0; ne < number_of_exponentials; ne++){
        double x = lifetime_spectrum[2 * ne];
        if(x == 0.0) continue;
        double lt_curr = lifetime_spectrum[2 * ne + 1];
        double tail_a = 1./(1.-exp(-period/lt_curr));
        double fit_curr = 0.;
        double exp_curr = std::exp(-dt/lt_curr);
        model[0] += dt_2 * instrument_response_function[0] * (exp_curr + 1.) * x;
        for(int i=convolution_start; i<convolution_stop; i++){
            fit_curr = (fit_curr + dt_2 * instrument_response_function[i - 1]) *
                       exp_curr + dt_2 * instrument_response_function[i];
            model[i] += fit_curr * x;
        }
        for(int i=convolution_stop; i<stop1; i++){
            fit_curr *= exp_curr;
            model[i] += fit_curr * x;
        }
        fit_curr *= exp(-(period_n - stop1) * dt / lt_curr);
        for(int i=0; i < convolution_stop; i++) {
            fit_curr *= exp_curr;
            model[i] += fit_curr * x * tail_a;
        }
    }
}


void fconv_cs_time_axis(
        double* output, int n_output,
        double* time_axis, int n_time_axis,
        double *instrument_response_function, int n_instrument_response_function,
        double* lifetime_spectrum, int n_lifetime_spectrum,
        int convolution_start,
        int convolution_stop
){
    int number_of_exponentials = n_lifetime_spectrum / 2;
    convolution_stop = convolution_stop > 0 ?
                       std::min(n_time_axis, std::min(n_instrument_response_function,
                                                      std::min(n_output, convolution_stop))) :
                       std::min(n_time_axis, std::min(n_instrument_response_function, n_output));
    convolution_start = std::max(convolution_start, 1);
#if VERBOSE
    std::clog << "convolve_lifetime_spectrum... " << std::endl;
    std::clog << "-- number_of_exponentials: " << number_of_exponentials << std::endl;
    std::clog << "-- convolution_start: " << convolution_start << std::endl;
    std::clog << "-- convolution_stop: " << convolution_stop << std::endl;
#endif
    for(int i=0; i<n_output; i++) output[i] = 0.0;
    for(int ne=0; ne<number_of_exponentials; ne++){
        double a = lifetime_spectrum[2 * ne];
        double current_lifetime = (lifetime_spectrum[2 * ne + 1]);
        if((a == 0.0) || (current_lifetime == 0.0)) continue;
        double current_model_value = 0.0;
        for(int i=convolution_start; i<convolution_stop; i++){
            double dt = dt = (time_axis[i] - time_axis[i - 1]);
            double dt_2 = dt / 2.0;
            double current_exponential = std::exp(-dt / current_lifetime);
            current_model_value = (current_model_value + dt_2 * instrument_response_function[i - 1]) *
                                  current_exponential + dt_2 * instrument_response_function[i];
            output[i] += current_model_value * a;
        }
    }
}

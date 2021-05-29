#include "include/DecayCurve.h"


void DecayCurve::add_arrays(
        double** output, int *n_output,
        double* curve1, int n_curve1,
        double* curve2, int n_curve2,
        double areal_fraction_curve2,
        int start,
        int stop
){
    *n_output = std::min(n_curve1, n_curve2);
    start = std::min(0, start);
    stop = stop < 0? *n_output: std::min(*n_output, stop);
    double sum_curve_1 = std::accumulate( curve1 + start, curve1 + stop, 1.0);
    double sum_curve_2 = std::accumulate( curve2 + start, curve2 + stop, 1.0);
    double f1 = (1. - areal_fraction_curve2);
    double f2 = areal_fraction_curve2 * sum_curve_1 / sum_curve_2;
    for(int i=start; i<stop;i++) curve1[i] *= f1;
    for(int i=start; i<stop;i++) curve1[i] += f2 * curve2[i];
    *output = curve1;
    *n_output = n_curve1;
}

std::vector<double> DecayCurve::shift_array(
        double* input, int n_input,
        double shift,
        bool set_outside,
        double outside_value
){
    auto out = std::vector<double>(n_input);
    if (n_input > 0) {
        // mod
        int r = (int) shift % n_input;
        int ts_i = r < 0 ? r + n_input : r;
        double ts_f = shift - std::floor(shift);
#if VERBOSE_FIT2X
        std::clog << "DecayCurve::shift_array" << std::endl;
        std::clog << "-- n_input: " << n_input << std::endl;
        std::clog << "-- shift: " << shift << std::endl;
        std::clog << "-- shift integer: " << ts_i << std::endl;
        std::clog << "-- shift float: " << ts_f << std::endl;
#endif
        std::vector<double> tmp1(input, input + n_input);
        std::vector<double> tmp2(input, input + n_input);
        std::rotate(tmp1.begin(), tmp1.begin() + ts_i, tmp1.end());
        std::rotate(tmp2.begin(), tmp2.begin() + (ts_i + 1), tmp2.end());
        for (int i = 0; i < n_input; i++) {
            out[i] = tmp1[i] * (1.0 - ts_f) + tmp2[i] * ts_f;
        }
        if (set_outside) {
            if (shift > 0) {
                for (int i = 0; i < r; i++)
                    out[i] = outside_value;
            } else if (shift < 0) {
                for (int i = n_input - r; i < n_input; i++)
                    out[i] = outside_value;
            }
        }

    }
    return out;
}

void DecayCurve::compute_noise(const char* noise_model){
    if(strcmp(noise_model, "poisson") == 0){
        for (size_t i = 0; i < size(); i++){
            ey[i] = std::max(1.0, std::sqrt(std::abs(y[i])));
        }
    }
}



DecayCurve& DecayCurve::operator+(const DecayCurve& other) const
{
    size_t n_max = std::min(size(), other.size());
    auto d = new DecayCurve();
    d->resize(n_max);
#if VERBOSE_FIT2X
    std::clog << "DecayCurve::operator+(const DecayCurve& other)" << std::endl;
    std::clog << "-- this->size():" << this->size() << std::endl;
    std::clog << "-- other->size():" << other.size() << std::endl;
    std::clog << "-- new->size():" << d->size() << std::endl;
#endif
    for(size_t i=0; i < n_max; i++){
        d->x[i] = x[i];
        d->y[i] = y[i] + other.y[i];
        d->ey[i] = std::sqrt(ey[i]*ey[i] + other.ey[i]*other.ey[i]);
    }
    d->acquisition_time = acquisition_time + other.acquisition_time;
    return *d;
}

DecayCurve& DecayCurve::operator*(const DecayCurve& other) const
{
    size_t n_max = std::min(size(), other.size());
    auto d = new DecayCurve();
    d->resize(n_max);
#if VERBOSE_FIT2X
    std::clog << "DecayCurve::operator*(const DecayCurve& other)" << std::endl;
    std::clog << "-- this->size():" << this->size() << std::endl;
    std::clog << "-- other->size():" << other.size() << std::endl;
    std::clog << "-- new->size():" << d->size() << std::endl;
#endif
    for(size_t i=0; i < n_max; i++){
        d->x[i] = x[i];
        d->y[i] = y[i] * other.y[i];
        double f1 = other.y[i]*ey[i];
        double f2 = y[i]*other.ey[i];
        d->ey[i] = std::sqrt(f1*f1 + f2*f2);
    }
    d->acquisition_time = acquisition_time + other.acquisition_time;
    return *d;
}

DecayCurve& DecayCurve::operator+(const double v) const {
    auto d = new DecayCurve();
    d->resize(size());
    for (size_t i = 0; i < size(); i++) {
        d->x[i] = x[i];
        d->y[i] = y[i] + v;
        d->ey[i] = ey[i];
    }
    return *d;
}

DecayCurve& DecayCurve::operator*(const double v) const {
    auto d = new DecayCurve();
    d->resize(size());
    for (size_t i = 0; i < size(); i++) {
        d->x[i] = x[i];
        d->y[i] = y[i] * v;
        d->ey[i] = ey[i] * v;
    }
    return *d;
}

DecayCurve& DecayCurve::operator+=(const double v) {
    for (size_t i = 0; i < size(); i++) {
        y[i] += v;
    }
    return *this;
}

DecayCurve& DecayCurve::operator*=(const double v) {
#if VERBOSE_FIT2X
    std::clog << "DecayCurve& DecayCurve::operator*=(const double v)" << std::endl;
    std::clog << "-- v:" << v << std::endl;
#endif
    for (size_t i = 0; i < size(); i++) {
        y[i] *= v;
        ey[i] *= v;
    }
    return *this;
}

DecayCurve& DecayCurve::operator=(const DecayCurve& other){
    // Guard self assignment
    if(this == &other){
        return *this;
    }

    // Assign other to this
    x = other.x;
    y = other.y; _y = other._y;
    ey = other.ey;

    return *this;
}

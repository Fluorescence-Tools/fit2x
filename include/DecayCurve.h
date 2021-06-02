#ifndef FIT2X_DECAYCURVE_H
#define FIT2X_DECAYCURVE_H

#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream> /* std::cerr */
#include <cstring> /* strcmp */
#include <memory> /* shared_ptr */
#include <limits> /* std::numeric_limits */

#include "tttrlib/TTTR.h"


class DecayCurve {

    friend class DecayConvolution;
    friend class DecayPileup;
    friend class DecayLinearization;
    friend class DecayScale;
    friend class DecayScore;
    friend class DecayPattern;

private:

    void compute_noise(const char* noise_model = "poisson");
    std::vector<double> _y; // original y values
    std::string noise_model = "NA";
    double current_shift = 0.0;
    double acquisition_time = 1e9;
    double avg_dx = 1.0;


    std::vector<double> x = std::vector<double>();
    std::vector<double> y = std::vector<double>();
    std::vector<double> ey = std::vector<double>();


public:

    static std::vector<double> shift_array(
            double *input, int n_input,
            double shift,
            bool set_outside = true,
            double outside_value = 0.0
    );

    /*!
     * Computes the sum of two arrays considering their respective
     * areal fraction.
     *
     * A weighted sum of two curves is computed. The weighted sum is
     * computed by the area of the first curve and the areal fraction
     * of the second curve. The area of the computed curve equals to
     * the area of the first input curve while the areal fraction of
     * the second input curve will be equal to the specified value in
     * the range specified by the input parameters.
     * modifies curve1 inplace
     *
     * @param output the computed output curve (array)
     * @param n_output the number of points in the output
     * @param curve1[in, out] the first input curve / array
     * @param n_curve1[in, out] number of points in the first array
     * @param curve2[in] second curve / array
     * @param n_curve2[in] number of points in the second array
     * @param areal_fraction_curve2[in] areal fraction of the second curve in the
     * output array
     * @param start[in] start index used for the area calculation
     * @param stop[in] stop index used for the area calculation
     */
    static void add_arrays(
            double **output_view, int *n_output,
            double *curve1, int n_curve1,
            double *curve2, int n_curve2,
            double areal_fraction_curve2,
            int start = 0,
            int stop = -1
    );

    size_t size() const{
        return x.size();
    }

    bool empty(){
        return x.empty();
    }

    std::vector<double> get_dx(){
        std::vector<double> dx(size(), 0.0);
        if(!empty()){
            for(size_t i = 1; i < x.size(); i++){
                dx[i] = x[i] - x[i - 1];
            }
        }
        return dx;
    }

    void resize(size_t n, double v=0.0, double dx=1.0){
        size_t old_size = size();
        // get last dx to extend the axis linearly
        auto dx_vec = get_dx();
        if(!dx_vec.empty()) dx = dx_vec[dx_vec.size() - 1];

        x.resize(n);
        y.resize(n, v); _y.resize(n, v);
        ey.resize(n, std::numeric_limits<double>::max());

        for(size_t i = old_size; i < n; i++){
            if(i > 0) x[i] = x[i - 1] + dx;
        }
    }


    double get_average_dx(){
        return avg_dx;
    }

    /*-------*/
    /* x     */
    /*-------*/
    std::vector<double>* get_x(){
        return &x;
    }

    void set_x(std::vector<double>& v){
        resize(v.size());
        x = v;
        // compute average dx
        auto dx = get_dx();
        avg_dx = std::accumulate(dx.begin(), dx.end(), 0.0) / (double) dx.size();
    }

    void set_x(double *input, int n_input){
        auto vec = std::vector<double>(input, input+n_input);
        set_x(vec);
    }

    /*-------*/
    /* y     */
    /*-------*/
    std::vector<double>* get_y(){
        return &y;
    }

    void set_y(std::vector<double>& v){
        resize(v.size());
        y = v; _y = v;
        compute_noise(noise_model.c_str());
    }

    void set_y(double *input, int n_input){
        auto vec = std::vector<double>(input, input+n_input);
        set_y(vec);
    }

    /*-------*/
    /* ey    */
    /*-------*/
    std::vector<double>* get_ey(){
        return &ey;
    }

    void set_ey(std::vector<double>& v){
        resize(v.size());
        ey = v;
    }

    void set_ey(double *input, int n_input){
        auto vec = std::vector<double>(input, input+n_input);
        set_ey(vec);
    }

    /*--------------------*/
    /* aquisition time    */
    /*--------------------*/
    void set_acquisition_time(double v) {
        if(acquisition_time < 0)
            acquisition_time = std::numeric_limits<double>::max();
        acquisition_time = v;
    }

    double get_acquisition_time() const {
        return acquisition_time;
    }

    /*--------------------*/
    /* time shift         */
    /*--------------------*/
    void set_shift(double v, double outside=0.0, bool shift_set_outside=true){
        current_shift = v;
        double* d = _y.data(); int n = _y.size();
        y = shift_array(d, n, current_shift, shift_set_outside, outside);
    }

    double get_shift(){
        return current_shift;
    }

    /*--------------------*/
    /* TTTR data          */
    /*--------------------*/
    void set_tttr(
            std::shared_ptr<TTTR> tttr_data,
            int tttr_micro_time_coarsening=1
    ){
        if(tttr_data != nullptr){
            double *hist; int n_hist;
            double *time; int n_time;
            TTTR::compute_microtime_histogram(
                    tttr_data.get(),
                    &hist, &n_hist,
                    &time, &n_time,
                    tttr_micro_time_coarsening
            );
            y.assign(hist, hist + n_hist);
            x.assign(time, time + n_time);
            free(hist); free(time);
            auto header = tttr_data->get_header();
            double micro_time_resolution =
                    header->get_micro_time_resolution() / tttr_micro_time_coarsening;
            for(size_t i = 0; i < x.size(); i++) x[i] = i * micro_time_resolution;
        }
    }

    DecayCurve(
            std::vector<double> x = std::vector<double>(),
            std::vector<double> y  = std::vector<double>(),
            std::vector<double> ey  = std::vector<double>(),
            double acquisition_time = -1,
            const char* noise_model = "poisson",
            int size = -1
            ){
        this->acquisition_time = acquisition_time;
        this->noise_model = noise_model;
        set_ey(ey);
        set_y(y);
        set_x(x);
        if(size > 0) resize(size);
    }

    DecayCurve& operator+(const DecayCurve& other) const;
    DecayCurve& operator*(const DecayCurve& other) const;
    DecayCurve& operator+(double v) const;
    DecayCurve& operator*(double v) const;
    DecayCurve& operator+=(double v);
    DecayCurve& operator*=(double v);
    DecayCurve& operator=(const DecayCurve& other);

};

#endif //FIT2X_DECAYCURVE_H

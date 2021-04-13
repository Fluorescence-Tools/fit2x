
#ifndef FIT2X_DECAYSCORE_H
#define FIT2X_DECAYSCORE_H

#include <string>
#include <iostream>

#include "DecayCurve.h"
#include "statistics.h"


class DecayScore
{

private:

    DecayCurve* _model = nullptr;
    DecayCurve* _data = nullptr;

    std::vector<double> _weighted_residuals;

    /// The range on which the scoring function is evaluated
    int _score_range_start = 0;
    int _score_range_stop = -1;

    int clip_by_size(int v, bool lower){
        int r = lower? 0: (int) _data->size();
        if(v > 0){
            r = std::min((int)_data->size(), v);
        }
        return r;
    }

    /// The used score type
    /*!
     * Options:
     */
    std::string _score_type= "poisson";

    void update_weighted_residuals() {
        int nm = _model->size();
        int nd = _data->size();
        int nw = std::min(nm, nd);
        _weighted_residuals.resize(nw);
#if VERBOSE_FIT2X
        std::clog << "Compute weighted residuals..." << std::endl;
        std::clog << "-- points in model function: " << nm << std::endl;
        std::clog << "-- points in data: " << nd << std::endl;
#endif
        for (int i = 0; i < nw; i++)
            _weighted_residuals[i] = (_data->y[i] - _model->y[i]) / _data->ey[i];
    }

public:

    void set_score_range_start(int v){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayScore::set_score_range_start" << std::endl;
#endif
        _score_range_start = v;
    }

    int get_score_range_start(){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayScore::get_score_range_start" << std::endl;
#endif
        return clip_by_size(_score_range_start, true);
    }

    void set_score_range_stop(int v){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayScore::set_score_range_stop" << std::endl;
#endif
        _score_range_stop = v;
    }

    int get_score_range_stop(){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayScore::get_score_range_stop" << std::endl;
#endif
        return clip_by_size(_score_range_stop, false);
    }

    void set_score_range(std::vector<int> v){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayScore::set_score_range" << std::endl;
#endif
        set_score_range_start(v[0]);
        set_score_range_stop(v[1]);
    }

    std::vector<int> get_score_range(){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayScore::DecayScore" << std::endl;
#endif
        return std::vector<int>({_score_range_start, _score_range_stop});
    }

    void set_model(DecayCurve* v){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayScore::set_model" << std::endl;
#endif
        _model = v;
    }

    void set_data(DecayCurve* v){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayScore::set_data" << std::endl;
#endif
        _data = v;
    }

    void set_score_type(std::string v){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayScore::set_score_type" << std::endl;
#endif
        _score_type = v;
    }

    std::string get_score_type(){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayScore::get_score_type" << std::endl;
#endif
        return _score_type;
    }

    void set(
            DecayCurve* model,
            DecayCurve* data,
            std::vector<int> score_range = std::vector<int>({0, -1}),
            std::string score_type = "poisson"
    ){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayScore::DecayScore" << std::endl;
#endif
        set_score_range_start(score_range[0]);
        set_score_range_stop(score_range[1]);
        set_score_type(score_type);
        set_model(model);
        set_data(data);
    }

    void get_weighted_residuals(double **output_view, int *n_output) {
#ifdef VERBOSE_FIT2X
        std::clog << "DecayScore::get_weighted_residuals" << std::endl;
#endif
        update_weighted_residuals();
        *n_output = _weighted_residuals.size();
        *output_view = _weighted_residuals.data();
    }

    double get_score(int x_min, int x_max, const char* score_type){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayScore::get_score" << std::endl;
#endif
#if VERBOSE_FIT2X
        std::clog << "CHI2" << std::endl;
#endif
        double v = 0.0;
        if((_model == nullptr) || (_data == nullptr)){
#if VERBOSE_FIT2X
            std::clog << "-- WARNING: decay or data not defined." << std::endl;
#endif
            v += INFINITY;
        } else {
            x_min = (x_min < 0) ?
                    get_score_range_start() :
                    clip_by_size(x_min, true);
            x_max = (x_max < 0) ?
                    get_score_range_stop():
                    clip_by_size(x_max, true);
#if VERBOSE_FIT2X
            std::clog << "-- data range: " << x_min << ", " << x_max << std::endl;
            std::clog << "-- score_type: " << score_type << std::endl;
#endif
            v = statistics::chi2_counting(
                    _data->y,
                    _model->y,
                    _data->w,
                    x_min, x_max,
                    score_type);
        }
        return v;
    }

    /// Evaluate and return the score
    double evaluate(){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayScore::evaluate" << std::endl;
#endif
        return get_score(
                get_score_range_start(), get_score_range_stop(),
                get_score_type().c_str()
        );
    }

    DecayScore(
            DecayCurve* model,
            DecayCurve* data,
            std::vector<int> score_range = std::vector<int>({0, -1}),
            std::string score_type= "poisson"
    ){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayScore::DecayScore" << std::endl;
#endif
        set(
                model,
                data,
                score_range,
                score_type
        );
    }

};

#endif //FIT2X_DECAYSCORE_H

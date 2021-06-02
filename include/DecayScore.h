
#ifndef FIT2X_DECAYSCORE_H
#define FIT2X_DECAYSCORE_H

#include <string>
#include <iostream>

#include "statistics.h"
#include "DecayRange.h"
#include "DecayCurve.h"


class DecayScore : public DecayRange {

private:

    DecayCurve* _model = nullptr;
    DecayCurve* _data = nullptr;
    DecayCurve* _default_model = nullptr;
    DecayCurve* _default_data = nullptr;
    std::vector<double> _weighted_residuals;
    std::string _score_type= "default";

public:

    DecayCurve* get_model(){
        auto re = _model;
        if(_model == nullptr){
            re = _default_model;
        }
        return re;
    }

    void set_model(DecayCurve* v){
        if(v != nullptr){
            _model = v;
        }
        get_data()->resize(get_model()->size());
    }

    DecayCurve* get_data(){
        auto re = _data;
        if(_data == nullptr){
            re = _default_data;
        }
        return re;
    }

    void set_data(DecayCurve* v){
        if(v != nullptr){
            _data = v;
        }
        get_model()->resize(get_data()->size());
    }

private:

    void update_weighted_residuals() {
        auto d = get_data();
        auto m = get_model();
        auto w = &_weighted_residuals;

        size_t start, stop;
        stop = get_stop(d);
        start = get_start(d);
        w->resize(stop - start);
        size_t j = 0;
        for (size_t i = start; i < stop; i++){
            (*w)[j] = (d->y[i] - m->y[i]) / d->ey[i];
            j++;
        }
#if VERBOSE_FIT2X
        std::clog << "Computed weighted residuals..." << std::endl;
        std::clog << "-- points in model function: " << m->size() << std::endl;
        std::clog << "-- points in data: " << d->size() << std::endl;
        std::clog << "-- weighted residual (w.r.) size: " << w->size() << std::endl;
        std::clog << "-- w.r. range: " << start << ", " << stop << std::endl;
#endif
    }

public:

    std::vector<double>& get_weighted_residuals() {
        update_weighted_residuals();
        return _weighted_residuals;
    }

    double get_score(int start=0, int stop=-1, const char* score_type = nullptr){
        // model data
        auto m = get_model();
        auto d = get_data();
        auto e = d->get_ey();
        // score range: start, stop
        start = std::min(get_start(d), get_start(m));
        stop = std::min(get_stop(d), get_stop(m));
        std::string st = (score_type == nullptr) ? get_score_type() : score_type;
        if(m->empty() || d->empty()){
            return INFINITY;
        }
#if VERBOSE_FIT2X
        std::clog << "DecayScore::get_score" << std::endl;
        std::clog << "-- score_type: " << st << std::endl;
        std::clog << "-- score range: " << start << ", " << stop << std::endl;
        std::clog << "-- model.size(): " << m->size() << std::endl;
        std::clog << "-- data.size(): " << d->size() << std::endl;
#endif
        double v = statistics::chi2_counting(
                *(d->get_y()), *(m->get_y()), *(e),
                get_start(d), get_stop(d),
                st.c_str()
        );
#if VERBOSE_FIT2X
        std::clog << "-- score: " << v << std::endl;
#endif
        return v;
    }

    void set_score_type(std::string v){
#if VERBOSE_FIT2X
        std::clog << "DecayScore::set_score_type" << std::endl;
#endif
        _score_type = v;
    }

    std::string get_score_type(){
#if VERBOSE_FIT2X
        std::clog << "DecayScore::get_score_type" << std::endl;
#endif
        return _score_type;
    }

    void set(
            DecayCurve* model, DecayCurve* data,
            std::string score_type = "poisson",
            int start=0, int stop=-1
    ){
#if VERBOSE_FIT2X
        std::clog << "DecayScore::DecayScore" << std::endl;
#endif
        set_start(start);
        set_stop(stop);
        set_score_type(score_type);
        set_model(model);
        set_data(data);
    }


    /// Evaluate and return the score
    double score(){
#if VERBOSE_FIT2X
        std::clog << "DecayScore::evaluate" << std::endl;
#endif
        return get_score(
                get_start(get_data()), get_stop(get_data()),
                get_score_type().c_str()
        );
    }

    DecayScore(
            DecayCurve* model= nullptr,
            DecayCurve* data= nullptr,
            std::string score_type= "default",
            int start=0, int stop=-1
    ) : DecayRange(start, stop){
#if VERBOSE_FIT2X
        std::clog << "DecayScore::DecayScore" << std::endl;
#endif
        _default_model = new DecayCurve();
        _default_data = new DecayCurve();
        _weighted_residuals = std::vector<double>();
        set(model, data, score_type);
    }

     ~DecayScore() override {
        delete _default_data;
        delete _default_model;
    }

};

#endif //FIT2X_DECAYSCORE_H

#ifndef FIT2X_DECAYPILEUP_H
#define FIT2X_DECAYPILEUP_H

#include <memory> /* shared_ptr */
#include <string>
#include <iostream> /* std::cerr */

#include "DecayCurve.h"
#include "DecayModifier.h"
#include "fsconv.h" /* add_pile_up_to_model */


class DecayPileup: public DecayModifier {

    friend class DecayCurve;

private:

    /// Dead time of the instrument in units of the lifetime (usually nanoseconds)
    double instrument_dead_time = std::numeric_limits<double>::epsilon();

    double repetition_rate = 100.0;

    /// If true pile up distortion will added to model function
    bool use_pile_up_correction = false;

    std::string pile_up_model = "coates";

    DecayCurve* data;

public:

    DecayCurve* get_data(){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayPileup::get_data" << std::endl;
#endif

        return data;
    }

    void set_data(DecayCurve* d){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayPileup::set_data" << std::endl;
#endif
        data = d;
    }

    void set_pile_up_model(std::string v){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayPileup::set_pile_up_model" << std::endl;
#endif
        pile_up_model = v;
    }

    std::string get_pile_up_model(){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayPileup::set_instrument_dead_time" << std::endl;
#endif
        return pile_up_model;
    }

    void set_instrument_dead_time(double v){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayPileup::set_instrument_dead_time" << std::endl;
#endif
        instrument_dead_time = v;
    }

    double get_instrument_dead_time(){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayPileup::get_instrument_dead_time" << std::endl;
#endif
        return instrument_dead_time;
    }

    void set_use_pile_up_correction(bool v){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayPileup::set_use_pile_up_correction" << std::endl;
#endif
        use_pile_up_correction = v;
    }

    bool get_use_pile_up_correction(){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayPileup::get_use_pile_up_correction" << std::endl;
#endif
        return use_pile_up_correction;
    }

    double get_repetition_rate(){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayPileup::get_repetition_rate" << std::endl;
#endif
        return repetition_rate;
    }

    void set_repetition_rate(double v){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayPileup::set_repetition_rate" << std::endl;
#endif
        repetition_rate = v;
    }

    DecayPileup(
            DecayCurve* data,
            const char* pile_up_model = "coates",
            double repetition_rate = 100,
            double instrument_dead_time = 120,
            bool use_pile_up_correction = true
    ){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayPileup::DecayPileup" << std::endl;
#endif

        set_data(data);
        set_pile_up_model(pile_up_model);
        set_repetition_rate(repetition_rate);
        set_instrument_dead_time(instrument_dead_time);
        set_use_pile_up_correction(use_pile_up_correction);
    }

    ~DecayPileup() = default;

    void resize(size_t n){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayPileup::resize" << std::endl;
#endif
        // For potential future use
    }

    void add(DecayCurve* model){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayPileup::add" << std::endl;
#endif
        if(use_pile_up_correction){
            if(data->equal_size(*model)){
                add_pile_up_to_model(
                        model->y.data(), model->y.size(),
                        data->y.data(), data->y.size(),
                        get_repetition_rate(),
                        get_instrument_dead_time(),
                        data->get_acquisition_time(),
                        pile_up_model.c_str()
                );
            } else{
                std::cerr << "WARNING: Data and model "
                             "differ in size. Pileup "
                             "not applied" << std::endl;
            }
        }
    }

};


#endif //FIT2X_DECAYPILEUP_H

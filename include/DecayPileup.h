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
        return data;
    }

    void set_data(DecayCurve* d){
        data = d;
    }

    void set_pile_up_model(std::string v){
        pile_up_model = v;
    }

    std::string get_pile_up_model(){
        return pile_up_model;
    }

    void set_instrument_dead_time(double v){
        instrument_dead_time = v;
    }

    double get_instrument_dead_time(){
        return instrument_dead_time;
    }

    void set_use_pile_up_correction(bool v){
        use_pile_up_correction = v;
    }

    bool get_use_pile_up_correction(){
        return use_pile_up_correction;
    }

    double get_repetition_rate(){
        return repetition_rate;
    }

    void set_repetition_rate(double v){
        repetition_rate = v;
    }

    DecayPileup(
            DecayCurve* data,
            const char* pile_up_model = "coates",
            double repetition_rate = 100,
            double instrument_dead_time = 120,
            bool use_pile_up_correction = true
    ){
        set_data(data);
        set_pile_up_model(pile_up_model);
        set_repetition_rate(repetition_rate);
        set_instrument_dead_time(instrument_dead_time);
        set_use_pile_up_correction(use_pile_up_correction);
    }

    ~DecayPileup() = default;

    void resize(size_t n){
        // For potential future use
    }

    void add(DecayCurve* model){
#if VERBOSE_FIT2X
        std::clog << "DecayPileup::add" << std::endl;
        std::clog << "-- use_pile_up_correction: " << use_pile_up_correction << std::endl;
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

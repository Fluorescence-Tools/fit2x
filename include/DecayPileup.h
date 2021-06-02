#ifndef FIT2X_DECAYPILEUP_H
#define FIT2X_DECAYPILEUP_H

#include <memory> /* shared_ptr */
#include <string>
#include <iostream> /* std::cerr */

#include "DecayCurve.h"
#include "DecayModifier.h"
#include "fsconv.h" /* add_pile_up_to_model */


class DecayPileup : public DecayModifier {

private:

    /// Dead time of the instrument in units of the lifetime (usually nanoseconds)
    double instrument_dead_time = std::numeric_limits<double>::epsilon();

    /// Repetition rate of the light source in MHz
    double repetition_rate = 100.0;

    /// Identifier of pile up model
    std::string pile_up_model = "coates";

public:

    void set_pile_up_model(std::string v){
        pile_up_model = v;
    }

    std::string get_pile_up_model(){
        return pile_up_model;
    }

    void set_repetition_rate(double v){
        repetition_rate = v;
    }

    double get_repetition_rate(){
        return repetition_rate;
    }

    void set_instrument_dead_time(double v){
        instrument_dead_time = std::abs(v);
    }

    double get_instrument_dead_time(){
        return instrument_dead_time;
    }

public:


    void add(DecayCurve* model){
        if(is_active()){
            model->resize(data->size());
            add_pile_up_to_model(
                    model->y.data(), model->y.size(),
                    data->y.data(), data->y.size(),
                    get_repetition_rate(),
                    get_instrument_dead_time(),
                    data->get_acquisition_time(),
                    pile_up_model.c_str()
            );
        }
    }

    DecayPileup(
        DecayCurve* data = nullptr,
        const char* pile_up_model = "coates",
        double repetition_rate = 100,
        double instrument_dead_time = 120,
        int start = 0,
        int stop = -1,
        bool active = true
    ) : DecayModifier(data, start, stop, active) {
#if VERBOSE_FIT2X
        std::clog << "DecayPileup::DecayPileup" << std::endl;
#endif
        set_pile_up_model(pile_up_model);
        set_repetition_rate(repetition_rate);
        set_instrument_dead_time(instrument_dead_time);
    }

};


#endif //FIT2X_DECAYPILEUP_H

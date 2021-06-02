#ifndef FIT2X_DecayPattern_H
#define FIT2X_DecayPattern_H

#include <iostream>
#include <algorithm>

# include "DecayModifier.h"

class DecayPattern: public DecayModifier{

private:

    /// The constant offset in the model
    double constant_offset = 0.0;

    /// Fraction (area) of background pattern
    double pattern_fraction = 0.0;

public:

    double get_pattern_fraction() const{
        return pattern_fraction;
    }

    void set_pattern_fraction(double v){
        pattern_fraction = std::max(0.0, std::min(1.0, std::abs(v)));
    }

    void set_pattern(DecayCurve* v){
        set_data(v);
    }

    DecayCurve* get_pattern(){
        return get_data();
    }

    void set_constant_offset(double v){
        constant_offset = v;
    }

    double get_constant_offset() const{
        return constant_offset;
    }

    void add(DecayCurve* decay) override{
        if(is_active()) {
            // resize background pattern to size of input decay
            resize(decay->size());

            int start = get_start(decay);
            int stop = get_stop(decay);

            double f = get_pattern_fraction();
            auto bg = get_pattern();


            // Compute pre-factors of decay and added pattern
            double ds, bs;
            if (f > 0.0) {
                ds = std::accumulate(decay->y.begin() + start, decay->y.begin() + stop, 0.0);
                bs = std::accumulate(bg->y.begin() + start, bg->y.begin() + stop, 0.0);
            } else {
                ds = bs = 1.0;
            }

            double fd = (1. - f);
            double fb = f * ds / bs;

            // mix constant offset, decay, and pattern
            double offset = get_constant_offset();
            for (int i = start; i < stop; i++) {
                decay->y[i] = decay->y[i] * fd + bg->y[i] * fb + offset;
            }
        }
    }

    DecayPattern(
            double constant_offset = 0.0,
            DecayCurve*pattern = nullptr,
            double pattern_fraction = 0.0,
            int start = 0,
            int stop = -1,
            bool active = true
    ) : DecayModifier(pattern, start, stop, active) {
#if VERBOSE_FIT2X
        std::clog << "DecayPattern::DecayPattern" << std::endl;
#endif
        set_pattern_fraction(pattern_fraction);
        set_constant_offset(constant_offset);
    }

};

#endif //FIT2X_DecayPattern_H
